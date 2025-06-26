using Quasar: QasmExpression, head, ClassicalVariable, GateDefinition, FunctionDefinition

# Define a type to represent measurement return values for each simulation path
struct MeasurementReturn
    path_measurements::Dict{Int, Vector{Qubit}}
end

###############################################
# Helper functions for evaluating expressions #
###############################################

"""
	evaluate_qubits(sim::BranchedSimulator, qubit_expr::QasmExpression) -> Vector{Vector{Qubit}}

Evaluate qubit expressions to get Qubit objects. This replaces the visitor's evaluate_qubits function.
Returns a list of list of Qubit objects. Each row represents all qubits to be used for a given active path.
The Qubit objects contain information about whether they've been measured.
"""
function evaluate_qubits(sim::BranchedSimulator, qubit_expr::QasmExpression)
	expr_type = head(qubit_expr)

	# Only process for active paths
	if isempty(sim.active_paths)
		return Vector{Vector{Qubit}}()
	end

	# Process all active paths
	all_qubits = Vector{Vector{Qubit}}()

	# Calculates the necessary qubits to be evaluated for each path and then appends the
	# result to the current list.
	for path_idx in sim.active_paths
		push!(all_qubits, evaluate_qubits_single_path(sim, path_idx, qubit_expr))
	end

	return all_qubits
end


function evaluate_qubits_single_path(sim::BranchedSimulator, path_idx::Int, qubit_expr::QasmExpression)
	all_qubits = Qubit[]
	expr_type = head(qubit_expr)

	# Ensure qubit_mapping exists for this path
	if path_idx > length(sim.qubit_mapping)
		push!(sim.qubit_mapping, Dict{String, Union{Qubit, Vector{Qubit}}}())
	end

	if expr_type == :identifier
		qubit_name = Quasar.name(qubit_expr)
		haskey(sim.qubit_mapping[path_idx], qubit_name) || error("Missing qubit '$qubit_name' for path $path_idx")

		value = sim.qubit_mapping[path_idx][qubit_name]

		if value isa Qubit # Single qubit
			push!(all_qubits, value)
		elseif value isa Vector{Qubit} # Qubit register
			append!(all_qubits, value)
		end
	elseif expr_type == :indexed_identifier
		qubit_name = Quasar.name(qubit_expr)

		# Get the index
		qubit_ix = _evolve_branched_ast(sim, qubit_expr.args[2])

		# Check if the register exists directly
		if haskey(sim.qubit_mapping[path_idx], qubit_name) &&
		   sim.qubit_mapping[path_idx][qubit_name] isa Vector{Qubit}
			# Get the register
			register = sim.qubit_mapping[path_idx][qubit_name]

			# Handle array indexing
			if qubit_ix isa Vector
				for ix in qubit_ix
					# Use 1-based indexing for register access
					if ix <= length(register)
						push!(all_qubits, register[ix])
					else
						error("Index $ix out of bounds for register '$qubit_name' of size $(length(register))")
					end
				end
			else
				# Use 1-based indexing for register access
				if qubit_ix <= length(register)
					push!(all_qubits, register[qubit_ix])
				else
					error("Index $qubit_ix out of bounds for register '$qubit_name' of size $(length(register))")
				end
			end
		else
			# Fall back to the old string-based approach for backward compatibility
			if qubit_ix isa Vector
				for ix in qubit_ix
					indexed_name = "$qubit_name[$ix]"
					haskey(sim.qubit_mapping[path_idx], indexed_name) || error("Missing qubit '$indexed_name' for path $path_idx")
					value = sim.qubit_mapping[path_idx][indexed_name]
					if value isa Qubit
						push!(all_qubits, value)
					elseif value isa Vector{Qubit} && length(value) == 1
						push!(all_qubits, value[1])
					end
				end
			else
				indexed_name = "$qubit_name[$qubit_ix]"
				haskey(sim.qubit_mapping[path_idx], indexed_name) || error("Missing qubit '$indexed_name' for path $path_idx")
				value = sim.qubit_mapping[path_idx][indexed_name]
				if value isa Qubit
					push!(all_qubits, value)
				elseif value isa Vector{Qubit} && length(value) == 1
					push!(all_qubits, value[1])
				end
			end
		end
	elseif expr_type == :array_literal
		# For array literals, recursively evaluate each element
		for element in qubit_expr.args
			append!(all_qubits, evaluate_qubits_single_path(sim, path_idx, element))
		end
	elseif expr_type == :hw_qubit
		# Hardware qubit reference (e.g., $0, $1)
		qubit_idx = parse(Int, replace(qubit_expr.args[1], "\$" => ""))
		
		# Get the qubit object
		qubit_obj = get_qubit_object(sim, path_idx, qubit_idx)
		push!(all_qubits, qubit_obj)
	end

	return all_qubits
end


"""
	evaluate_binary_op(op::Symbol, lhs, rhs)

Evaluate binary operations. This is a direct copy of the visitor's function.
"""
function evaluate_binary_op(op::Symbol, lhs, rhs)
	op == :< && return lhs < rhs
	op == :> && return lhs > rhs
	op == :<= && return lhs <= rhs
	op == :>= && return lhs >= rhs
	op == Symbol("=") && return rhs
	op == Symbol("!=") && return lhs != rhs
	op == Symbol("==") && return lhs == rhs
	op == :+ && return lhs + rhs
	op == :- && return lhs - rhs
	op == :* && return lhs * rhs
	op == :/ && return lhs / rhs
	op == :% && return lhs % rhs
	op == Symbol("<<") && return lhs << rhs
	op == Symbol(">>") && return lhs >> rhs
	op == Symbol("**") && return lhs ^ rhs
	op == Symbol("&&") && return lhs && rhs
	op == Symbol("||") && return lhs || rhs
	op == :| && return lhs .| rhs
	op == :& && return lhs .& rhs
	op == :^ && return lhs .⊻ rhs
	error("Unknown binary operator: $op")
end

"""
	evaluate_unary_op(op::Symbol, arg)

Evaluate unary operations. This is a direct copy of the visitor's function.
"""
function evaluate_unary_op(op::Symbol, arg)
	op == :! && return !arg
	op == :~ && return .!arg
	op == :- && return -arg
	error("Unknown unary operator: $op")
end

"""
	evaluate_modifiers(sim::BranchedSimulator, expr::QasmExpression)

Evaluate gate modifiers. This replaces the visitor's evaluate_modifiers function.
"""
function evaluate_modifiers(sim::BranchedSimulator, expr::QasmExpression)
	if head(expr) == :power_mod
		pow_expr = QasmExpression(:pow, _evolve_branched_ast(sim, expr.args[1]))
		return (pow_expr, expr.args[2])
	elseif head(expr) == :inverse_mod
		return (QasmExpression(:inv), expr.args[1])
	elseif head(expr) ∈ (:control_mod, :negctrl_mod)
		has_argument = length(expr.args) > 1
		if has_argument
			arg_val = _evolve_branched_ast(sim, first(expr.args))
			isinteger(arg_val) || error("Cannot apply non-integer ($arg_val) number of controls or negcontrols.")
			true_inner = expr.args[2]
			inner = QasmExpression(head(expr), true_inner)
			while arg_val > 2
				inner = QasmExpression(head(expr), inner)
				arg_val -= 1
			end
		else
			inner = expr.args[1]
		end
		new_head = head(expr) == :control_mod ? :ctrl : :negctrl
		return (QasmExpression(new_head), inner)
	else
		error("Unknown modifier type: $(head(expr))")
	end
end

#####################################
# Main Function for State Evolution #
#####################################

"""
	evolve_branched!(simulator::AbstractSimulator, program::QasmExpression, inputs::Dict{String, <:Any}) -> BranchedSimulator

Evolve a quantum program using a branched approach that handles measurements and control flow.
"""
function evolve_branched(simulator::AbstractSimulator, program::QasmExpression, inputs::Dict{String, <:Any})
	# Create a branched simulator with integrated visitor functionality
	branched_sim = BranchedSimulator(simulator; inputs = inputs)

	# Process the AST
	_evolve_branched_ast(branched_sim, program)

	return branched_sim
end


_evolve_branched_ast(sim::BranchedSimulator, i::Number) = i
_evolve_branched_ast(sim::BranchedSimulator, i::String) = i
_evolve_branched_ast(sim::BranchedSimulator, i::BitVector) = i
_evolve_branched_ast(sim::BranchedSimulator, i::NTuple{N, <:Number}) where {N} = i
_evolve_branched_ast(sim::BranchedSimulator, i::Vector{<:Number}) = i
"""
	_evolve_branched_ast(sim::BranchedSimulator, expr::QasmExpression)

Process an AST node in the branched simulation model.
"""
function _evolve_branched_ast(sim::BranchedSimulator, expr::QasmExpression)
	println(expr)
	println(sim.states)
	println(sim.variables)
	println(sim.qubit_mapping)
	println(sim.measurements)
	expr_type = head(expr)

	# Program node - process each child node in sequence
	if expr_type == :program
		for child_expr in expr.args
			head(child_expr) == :end && return
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && return  # No active paths left
		end

		# Scope node - process each child node in sequence
	elseif expr_type == :scope
		for child_expr in expr.args
			child_head = head(child_expr)
			if child_head in (:end, :continue, :break)
				return
			end
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && return  # No active paths left
		end

		# Simple operations that don't affect quantum state
	elseif expr_type == :version
		# Version node - ignore
		return

	elseif expr_type == :reset
		# Reset operation - reset qubit to |0⟩ state
		targets_by_path = evaluate_qubits(sim, expr.args[1])
		for (path_idx, targets) in zip(sim.active_paths, targets_by_path)
			for qubit in targets
				# Check if this qubit was previously measured
				if qubit.measured
					# Expand the state to reincorporate the qubit
					expand_state!(sim, path_idx, qubit.index)
				end
				
				# Apply reset to the qubit's current index in the state
				qubit_idx = get_qubit_index(sim, path_idx, qubit)
				_apply_reset!(sim.states[path_idx], qubit_idx)
			end
		end

	elseif expr_type == :barrier || expr_type == :delay
		# Barrier/delay operation - no effect on simulation
		return

	elseif expr_type in [:stretch, :duration, :box, :output]
		# Unsupported operations
		error("Simulator doesn't support $expr_type operations")

		# Variable declarations and control flow
	elseif expr_type == :input
		# Input parameter declaration
		var_name = Quasar.name(expr)
		var_type = expr.args[1].args[1]
		haskey(sim.inputs, var_name) || error("Missing input variable '$var_name'.")
		for path_idx in sim.active_paths
			var = ClassicalVariable(var_name, var_type, sim.inputs[var_name], true)
			sim.variables[path_idx][var_name] = var
		end

	elseif expr_type ∈ (:continue, :break)
		# Loop control statements
		return

		# For loop
	elseif expr_type == :for
		_handle_for_loop(sim, expr)

		# Switch statement
	elseif expr_type == :switch
		_handle_switch_statement(sim, expr)

		# Alias declaration
	elseif expr_type == :alias
		alias_name = Quasar.name(expr)
		right_hand_side = expr.args[1].args[1].args[end]

		# Process the alias based on the right-hand side type
		if head(right_hand_side) == :binary_op
			# Handle concatenation
			right_hand_side.args[1] == Symbol("++") || error("Right hand side of alias must be either an identifier, indexed identifier, or concatenation")
			concat_left = right_hand_side.args[2]
			concat_right = right_hand_side.args[3]

			# Check if both sides are qubits
			is_left_qubit = false
			is_right_qubit = false

			# Only process for active paths
			for path_idx in sim.active_paths
				# Ensure qubit_mapping exists for this path
				if path_idx > length(sim.qubit_mapping)
					push!(sim.qubit_mapping, Dict{String, Union{Qubit, Vector{Qubit}}}())
				end

				left_name = Quasar.name(concat_left)
				right_name = Quasar.name(concat_right)

				is_left_qubit = haskey(sim.qubit_mapping[path_idx], left_name)
				is_right_qubit = haskey(sim.qubit_mapping[path_idx], right_name)

				# Check if trying to concatenate qubit and classical arrays
				(is_left_qubit ⊻ is_right_qubit) && error("Cannot concatenate qubit and classical arrays")

				if is_left_qubit
					# Check if either side is a physical qubit
					if head(concat_left) == :hw_qubit || head(concat_right) == :hw_qubit
						error("Cannot alias physical qubits. The let keyword allows declared quantum bits and registers to be referred to by another name.")
					end

					# Both are qubits, concatenate them
					left_qs_by_path = evaluate_qubits(sim, concat_left)
					right_qs_by_path = evaluate_qubits(sim, concat_right)
					
					# Since we're already in a loop for each path_idx, we can directly access the current path's qubits
					left_qs = left_qs_by_path[findfirst(p -> p == path_idx, sim.active_paths)]
					right_qs = right_qs_by_path[findfirst(p -> p == path_idx, sim.active_paths)]
					alias_qubits = vcat(left_qs, right_qs)

					# Store the concatenated register
					sim.qubit_mapping[path_idx][alias_name] = alias_qubits

					# For backward compatibility, also map individual qubits
					for (i, qubit) in enumerate(alias_qubits)
						sim.qubit_mapping[path_idx]["$alias_name[$(i-1)]"] = qubit
					end
				else
					# Both are classical variables, which is not allowed
					error("The let keyword only allows aliasing for quantum bits and registers, not classical variables.")
				end
			end
		elseif head(right_hand_side) == :identifier
			referent_name = Quasar.name(right_hand_side)

			# Only process for active paths
			for path_idx in sim.active_paths
				# Ensure qubit_mapping exists for this path
				if path_idx > length(sim.qubit_mapping)
					push!(sim.qubit_mapping, Dict{String, Union{Qubit, Vector{Qubit}}}())
				end

				# Check if it's a qubit alias
				if haskey(sim.qubit_mapping[path_idx], referent_name)
					# Check if it's a physical qubit
					if head(right_hand_side) == :hw_qubit
						error("Cannot alias physical qubits. The let keyword allows declared quantum bits and registers to be referred to by another name.")
					end

					# Copy the qubit mapping
					sim.qubit_mapping[path_idx][alias_name] = sim.qubit_mapping[path_idx][referent_name]

					# For backward compatibility, also map individual qubits if it's a register
					if sim.qubit_mapping[path_idx][referent_name] isa Vector{Qubit}
						qubit_register = sim.qubit_mapping[path_idx][referent_name]
						for (i, qubit) in enumerate(qubit_register)
							sim.qubit_mapping[path_idx]["$alias_name[$(i-1)]"] = qubit
						end
					end
				else
					# It's a classical variable alias, which is not allowed
					error("The let keyword only allows aliasing for quantum bits and registers, not classical variables.")
				end
			end
		elseif head(right_hand_side) == :indexed_identifier
			referent_name = Quasar.name(right_hand_side)

			# Only process for active paths
			for path_idx in sim.active_paths
				# Ensure qubit_mapping exists for this path
				if path_idx > length(sim.qubit_mapping)
					push!(sim.qubit_mapping, Dict{String, Union{Qubit, Vector{Qubit}}}())
				end

				# Check if it's a qubit alias
				if haskey(sim.qubit_mapping[path_idx], referent_name)
					# Get the qubits from the indexed identifier
					qubits_by_path = evaluate_qubits(sim, right_hand_side)
					# Since we're already in a loop for each path_idx, we can directly access the current path's qubits
					alias_qubits = qubits_by_path[findfirst(p -> p == path_idx, sim.active_paths)]

					# Store the register
					sim.qubit_mapping[path_idx][alias_name] = alias_qubits

					# For backward compatibility, also map individual qubits
					for (i, qubit) in enumerate(alias_qubits)
						sim.qubit_mapping[path_idx]["$alias_name[$(i-1)]"] = qubit
					end
				else
					# It's a classical variable alias, which is not allowed
					error("The let keyword only allows aliasing for quantum bits and registers, not classical variables.")
				end
			end
		elseif head(right_hand_side) == :hw_qubit
			# Trying to alias a physical qubit
			error("Cannot alias physical qubits. The let keyword allows declared quantum bits and registers to be referred to by another name.")
		else
			error("Right hand side of alias must be either an identifier, indexed identifier, or concatenation")
		end

		# Identifier reference
	elseif expr_type == :identifier
		id_name = Quasar.name(expr)

		# Only process for active paths
		if isempty(sim.active_paths)
			error("No active paths when evaluating identifier $id_name")
		end

		# Get the current path
		path_idx = first(sim.active_paths)

		# First check path-specific variables
		var_value = get_variable(sim, path_idx, id_name)
		if !isnothing(var_value)
			return var_value
		end

		# Then check path-specific qubit_mapping
		if path_idx <= length(sim.qubit_mapping) &&
		   haskey(sim.qubit_mapping[path_idx], id_name)
			return sim.qubit_mapping[path_idx][id_name]
		end

		# If we get here, the identifier is not defined
		error("No identifier $id_name defined for path $path_idx")

		# Indexed identifier reference
	elseif expr_type == :indexed_identifier
		identifier_name = Quasar.name(expr)

		# Only process for active paths
		if isempty(sim.active_paths)
			error("No active paths when evaluating indexed identifier $identifier_name")
		end

		# Get the current path
		path_idx = first(sim.active_paths)

		# Get the index
		ix = _evolve_branched_ast(sim, expr.args[2])
		flat_ix = (ix isa Vector) ? ix : [ix]
		flat_ix = flat_ix .+ 1  # Convert to 1-based indexing

		# First check path-specific variables
		var_value = get_variable(sim, path_idx, identifier_name)
		if !isnothing(var_value) && var_value isa Vector
			return var_value[only(flat_ix)]
		end

		# Then check path-specific qubit_mapping
		if path_idx <= length(sim.qubit_mapping)
			return evaluate_qubits(sim, expr)
		end

		# If we get here, the indexed identifier is not defined
		error("No indexed identifier $identifier_name defined for path $path_idx")

		# Conditional statement
	elseif expr_type == :if
		_handle_conditional(sim, expr)

		# While loop
	elseif expr_type == :while
		_handle_while_loop(sim, expr)

		# Classical variable assignment
	elseif expr_type == :classical_assignment
		_handle_classical_assignment(sim, expr)

		# Classical variable declaration
	elseif expr_type == :classical_declaration
		_handle_classical_declaration(sim, expr)

		# Constant declaration
	elseif expr_type == :const_declaration
		_handle_const_declaration(sim, expr)

		# Qubit declaration
	elseif expr_type == :qubit_declaration
		_handle_qubit_declaration(sim, expr)

		# Gate modifiers (power, inverse, control, negctrl)
	elseif expr_type ∈ [:power_mod, :inverse_mod, :control_mod, :negctrl_mod]
		_handle_gate_modifiers(sim, expr)

		# Gate call
	elseif expr_type == :gate_call
		_handle_gate_call(sim, expr)

		# Gate definition
	elseif expr_type == :gate_definition
		_handle_gate_definition(sim, expr)

		# Measurement operation - key part that creates path branching
	elseif expr_type == :measure
		return _handle_measurement(sim, expr)

		# Function call
	elseif expr_type == :function_call
		_handle_function_call(sim, expr)

		# Function definition
	elseif expr_type == :function_definition
		_handle_function_definition(sim, expr)

		# Pragma directive
	elseif expr_type == :pragma
		# For simplicity, ignore pragmas
		return

		# Literal values
	elseif expr_type ∈ [:integer_literal, :float_literal, :string_literal, :complex_literal, :irrational_literal, :boolean_literal, :duration_literal]
		# Return the literal value
		return expr.args[1]

		# Array literal
	elseif expr_type == :array_literal
		# Evaluate array elements
		return [_evolve_branched_ast(sim, element) for element in expr.args]

		# Range expression
	elseif expr_type == :range
		# Evaluate range parameters
		start = _evolve_branched_ast(sim, expr.args[1])
		stop = _evolve_branched_ast(sim, expr.args[2])
		step = length(expr.args) > 2 ? _evolve_branched_ast(sim, expr.args[3]) : 1
		return collect(start:step:stop)

		# Expression evaluation
	elseif expr_type == :binary_op
		# Binary operation
		op = expr.args[1]
		lhs = _evolve_branched_ast(sim, expr.args[2])
		rhs = _evolve_branched_ast(sim, expr.args[3])

		# Evaluate based on operator type
		return evaluate_binary_op(op, lhs, rhs)

	elseif expr_type == :unary_op
		# Unary operation
		op = expr.args[1]
		arg = _evolve_branched_ast(sim, expr.args[2])

		# Evaluate using existing function
		return evaluate_unary_op(op, arg)

		# Type cast
	elseif expr_type == :cast
		casting_to = expr.args[1].args[1]
		value = _evolve_branched_ast(sim, expr.args[2])

		# Handle basic type conversions
		if casting_to == Bool
			return value > 0
		elseif casting_to == Int
			return Int(value)
		elseif casting_to == Float64
			return Float64(value)
		else
			error("Unsupported cast type: $casting_to")
		end

		# Unsupported expression
	else
		error("Cannot process expression of type $expr_type.")
	end
end

################################
# Handlers for State Evolution #
################################

"""
	_handle_measurement(sim::BranchedSimulator, expr::QasmExpression)

Handle a measurement operation, creating branches for different outcomes.
Uses the memory optimization feature to reduce state size after measurement.
Returns a MeasurementReturn object containing the measured qubits for each path.
"""
function _handle_measurement(sim::BranchedSimulator, expr::QasmExpression)
	# TEST HERE
	sim.states[1] = [0.5+0im, 0.5+0im, 0.5+0im, 0.5+0im]

	# Get qubits to measure for each path
	qubits_by_path = evaluate_qubits(sim, expr.args[1])

	# Now measure all qubits for all paths
	new_active_paths = Int[]

	# Structure to store measured qubits for each path
	path_measurements = Dict{Int, Vector{Qubit}}()
	
	for (path_idx_id, qubits_to_measure) in enumerate(qubits_by_path)
		path_idx = sim.active_paths[path_idx_id]

		# Everytime a measurement occurs, the resulting path indices are added to the following list
		paths_to_measure = [path_idx]

		# Initialize path_measurements for the initial path
		path_measurements[path_idx] = Qubit[]

		for qubit in qubits_to_measure
			# Skip already measured qubits for actual measurement
			if qubit.measured
				# Still record the qubit in all paths
				for idx in paths_to_measure
					# Create a deep copy of the qubit to avoid aliasing issues
					push!(path_measurements[idx], Qubit(qubit.index, qubit.name, qubit.measured))
				end
				continue
			end
			
			# Debug: Print qubit info before measurement
			println("Measuring qubit: $(qubit.index), name: $(qubit.name), measured: $(qubit.measured)")
			
			# This inner for loop iterates through all sub paths created by a given path undergoing a measurement
			current_paths = copy(paths_to_measure)
			for idx in current_paths
				# Debug: Print state before measurement
				println("Path $idx state before measurement: $(sim.states[idx])")
				
				# Map the qubit to its current position in the state
				current_qubit_idx = get_qubit_index(sim, idx, qubit)
				println("Qubit $(qubit.index) mapped to position $current_qubit_idx in state")
				
				# Get measurement probabilities
				probs = get_measurement_probabilities(sim.states[idx], current_qubit_idx)
				println("Measurement probabilities: $probs")
				
				# Measure the qubit on this path
				path_split = measure_qubit(sim, idx, qubit)
				if idx == path_idx
					mark_unmeasured(qubit)
				end
				println("Path split: $path_split")
				
				# Get the updated qubit object after measurement
				measured_qubit = get_qubit_object(sim, idx, qubit.index)
				
				# Add the measured qubit to this path's measurements
				push!(path_measurements[idx], Qubit(measured_qubit.index, measured_qubit.name, measured_qubit.measured))
				
				if path_split
					# A new path was created
					added_idx = length(sim.states) # Since we add the index at the end of the states vector
					push!(paths_to_measure, added_idx)
					
					# Copy all previous measurements to the new path, creating new Qubit objects to avoid aliasing
					path_measurements[added_idx] = [Qubit(q.index, q.name, q.measured) for q in path_measurements[idx]]
					
					# Debug: Print new path info
					println("Created new path $added_idx with state: $(sim.states[added_idx])")
				end
			end
			mark_measured(qubit)
		end

		# Add all the new subpaths to the active paths list
		append!(new_active_paths, paths_to_measure)
	end
	
	# Update active paths all at once
	sim.active_paths = new_active_paths

	# Debug: Print final paths and states
	println("Final active paths: $(sim.active_paths)")
	for path_idx in sim.active_paths
		println("Path $path_idx state: $(sim.states[path_idx])")
		println("Path $path_idx measurements: $(sim.measurements[path_idx])")
	end

	# Return a MeasurementReturn object with the qubit measurements for each path
	return MeasurementReturn(path_measurements)
end

"""
	_handle_conditional(sim::BranchedSimulator, expr::QasmExpression)

Handle a conditional statement, filtering paths based on the condition.
"""
function _handle_conditional(sim::BranchedSimulator, expr::QasmExpression)
	# Find if there's an else branch
	has_else = findfirst(e->head(e) == :else, convert(Vector{QasmExpression}, expr.args))
	last_expr = !isnothing(has_else) ? length(expr.args) - 1 : length(expr.args)

	# Evaluate condition for each active path
	true_paths = Int[]
	false_paths = Int[]

	for path_idx in sim.active_paths
		# Evaluate condition in this path's context
		condition_value = _evaluate_condition(sim, path_idx, expr.args[1])

		if condition_value
			push!(true_paths, path_idx)
		else
			push!(false_paths, path_idx)
		end
	end

	# Process if branch for true paths
	if !isempty(true_paths)
		# Save original active paths
		original_active_paths = copy(sim.active_paths)

		# Set active paths to only those where condition is true
		sim.active_paths = true_paths

		# Process if branch
		for child_expr in expr.args[2:last_expr]
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && break  # No active paths left
		end

		# Restore active paths that were not processed in if branch
		append!(sim.active_paths, setdiff(original_active_paths, true_paths))
	end

	# Process else branch for false paths
	if !isnothing(has_else) && !isempty(false_paths)
		# Save original active paths
		original_active_paths = copy(sim.active_paths)

		# Set active paths to only those where condition is false
		sim.active_paths = false_paths

		# Process else branch
		for child_expr in expr.args[has_else].args
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && break  # No active paths left
		end

		# Restore active paths that were not processed in else branch
		append!(sim.active_paths, setdiff(original_active_paths, false_paths))
	end
end

"""
	_evaluate_condition(sim::BranchedSimulator, path_idx::Int, condition_expr::QasmExpression)

Evaluate a condition expression in the context of a specific path.
This handles path-specific variables and measurement results.
"""
function _evaluate_condition(sim::BranchedSimulator, path_idx::Int, condition_expr::QasmExpression)
	# Save the current active paths
	original_active_paths = copy(sim.active_paths)

	# Set active paths to only the current path
	sim.active_paths = [path_idx]

	# Handle different types of condition expressions
	if head(condition_expr) == :identifier
		# Simple variable reference
		var_name = name(condition_expr)

		# Check path-specific variables from the simulator
		result = get_variable(sim, path_idx, var_name)
		if !isnothing(result)
			# Restore active paths
			sim.active_paths = original_active_paths
			return result > 0
		end

		# If we get here, the variable is not defined
		# Restore active paths
		sim.active_paths = original_active_paths
		return false
	elseif head(condition_expr) == :binary_op
		# Binary operation (e.g., ==, !=, <, >, etc.)
		op = condition_expr.args[1]
		lhs = _evolve_branched_ast(sim, condition_expr.args[2])
		rhs = _evolve_branched_ast(sim, condition_expr.args[3])

		# Evaluate the operation
		result = evaluate_binary_op(op, lhs, rhs)

		# Restore active paths
		sim.active_paths = original_active_paths
		return result
	elseif head(condition_expr) == :indexed_identifier
		# Indexed variable reference (e.g., b[0])
		var_name = name(condition_expr)

		# Get the index
		ix = _evolve_branched_ast(sim, condition_expr.args[2])
		flat_ix = (ix isa Vector) ? ix : [ix]
		flat_ix = flat_ix .+ 1  # Convert to 1-based indexing

		# Check path-specific variables from the simulator
		var_value = get_variable(sim, path_idx, var_name)
		if !isnothing(var_value) && var_value isa Vector
			result = var_value[only(flat_ix)]
			# Restore active paths
			sim.active_paths = original_active_paths
			return result > 0
		end

		# If we get here, the variable is not defined or not a vector
		# Restore active paths
		sim.active_paths = original_active_paths
		return false
	else
		# For other expression types, evaluate using the AST processor
		result = _evolve_branched_ast(sim, condition_expr)

		# Restore active paths
		sim.active_paths = original_active_paths

		# Convert to boolean
		return result isa Bool ? result : result > 0
	end

	# Restore active paths
	sim.active_paths = original_active_paths

	# Default case
	return false
end

"""
	_handle_for_loop(sim::BranchedSimulator, expr::QasmExpression)

Handle a for loop in the branched simulation model.
"""
function _handle_for_loop(sim::BranchedSimulator, expr::QasmExpression)
	for_loop = convert(Vector{QasmExpression}, expr.args)
	loop_variable_type = for_loop[1].args[1]
	loop_variable_name = for_loop[2].args[1]::String

	# Evaluate loop range
	loop_variable_values = _evolve_branched_ast(sim, for_loop[3])
	loop_body = for_loop[4]::QasmExpression

	# For each value in the range
	for loop_value in loop_variable_values
		# Set the variable for each active path
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, loop_variable_name, ClassicalVariable(loop_variable_name, loop_variable_type, loop_value, false))
		end

		# Process loop body
		_evolve_branched_ast(sim, loop_body)

		# If no active paths left, break early
		isempty(sim.active_paths) && break
	end

	# Clean up from all paths
	for path_idx in sim.active_paths
		if haskey(sim.variables[path_idx], loop_variable_name)
			delete!(sim.variables[path_idx], loop_variable_name)
		end
	end
end

"""
	_handle_switch_statement(sim::BranchedSimulator, expr::QasmExpression)

Handle a switch statement in the branched simulation model.
"""
function _handle_switch_statement(sim::BranchedSimulator, expr::QasmExpression)
	all_cases = convert(Vector{QasmExpression}, expr.args[2:end])
	default_idx = findfirst(expr->head(expr) == :default, all_cases)

	# Group paths by their case values
	case_paths = Dict{Any, Vector{Int}}()
	default_paths = Int[]

	# For each active path, evaluate the condition and determine which case it should follow
	original_active_paths = copy(sim.active_paths)

	for path_idx in original_active_paths
		# Temporarily set this path as the only active one
		sim.active_paths = [path_idx]

		# Evaluate switch condition for this path
		case_val = _evolve_branched_ast(sim, expr.args[1])

		# Find matching case for this path
		path_matched = false

		for (case_idx, case) in enumerate(all_cases)
			if head(case) == :case
				# Evaluate case value for this path
				sim.active_paths = [path_idx]
				case_value_result = _evolve_branched_ast(sim, case.args[1])

				if case_val ∈ case_value_result
					# This path matches this case
					if !haskey(case_paths, case_idx)
						case_paths[case_idx] = Int[]
					end
					push!(case_paths[case_idx], path_idx)
					path_matched = true
					break
				end
			end
		end

		# If no case matched, add to default paths
		if !path_matched
			push!(default_paths, path_idx)
		end
	end

	# Save all active paths before processing cases
	all_active_paths = copy(sim.active_paths)

	# Collect paths that survive each case
	surviving_paths = Int[]

	# Process each case with its matching paths
	for (case_idx, paths) in case_paths
		# Skip if no paths match this case
		isempty(paths) && continue

		# Set active paths to only those that match this case
		sim.active_paths = paths

		# Process case body
		case = all_cases[case_idx]
		for child_expr in convert(Vector{QasmExpression}, case.args[2:end])
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && break  # No active paths left for this case
		end

		# Add surviving paths to the collection
		append!(surviving_paths, sim.active_paths)
	end

	# Process default case if needed
	if !isempty(default_paths) && !isnothing(default_idx)
		# Set active paths to only those that need the default case
		sim.active_paths = default_paths

		# Process default case body
		for child_expr in convert(Vector{QasmExpression}, all_cases[default_idx].args)
			_evolve_branched_ast(sim, child_expr)
			isempty(sim.active_paths) && break  # No active paths left for default
		end

		# Add surviving paths to the collection
		append!(surviving_paths, sim.active_paths)
	elseif !isempty(default_paths) && isnothing(default_idx)
		error("No case matched for some paths and no default defined.")
	end

	# Update active paths to include all paths that survived their respective cases
	sim.active_paths = surviving_paths

	# If no paths survived, restore original paths
	if isempty(sim.active_paths)
		sim.active_paths = all_active_paths
	end
end

"""
	_handle_while_loop(sim::BranchedSimulator, expr::QasmExpression)

Handle a while loop in the branched simulation model.
"""
function _handle_while_loop(sim::BranchedSimulator, expr::QasmExpression)
	loop_body = expr.args[2]

	# Keep track of paths that should continue looping
	continue_paths = copy(sim.active_paths)

	while !isempty(continue_paths)
		# Save original active paths
		original_active_paths = copy(sim.active_paths)

		# Set active paths to only those that should continue looping
		sim.active_paths = continue_paths

		# Evaluate condition for each path
		new_continue_paths = Int[]
		for path_idx in continue_paths
			# Temporarily set this path as the only active one
			sim.active_paths = [path_idx]

			# Evaluate condition
			condition_value = _evolve_branched_ast(sim, expr.args[1]) > 0

			if condition_value
				push!(new_continue_paths, path_idx)
			end
		end

		# If no paths should continue, break
		isempty(new_continue_paths) && break

		# Set active paths to only those that should continue
		sim.active_paths = new_continue_paths

		# Process loop body
		_evolve_branched_ast(sim, loop_body)

		# Update continue_paths for next iteration
		continue_paths = copy(sim.active_paths)

		# Restore paths that didn't enter the loop
		append!(sim.active_paths, setdiff(original_active_paths, new_continue_paths))
	end
end

"""
	_handle_classical_assignment(sim::BranchedSimulator, expr::QasmExpression)

Handle a classical variable assignment in the branched simulation model.
"""
function _handle_classical_assignment(sim::BranchedSimulator, expr::QasmExpression)
	op = expr.args[1].args[1]
	lhs = expr.args[1].args[2]
	rhs = expr.args[1].args[3]

	# Evaluate the right-hand side first
	# This will handle any measurements or complex expressions
	rhs_value = _evolve_branched_ast(sim, rhs)

	var_name = Quasar.name(lhs)
	
	# Check if the right-hand side is a measurement result
	is_measurement = rhs_value isa MeasurementReturn

	# This is for defining a register
	if head(lhs) == :identifier
		if is_measurement
			# For each measured path, assign the corresponding classical register bit values
			for (path_idx_index, path_idx) in enumerate(sim.active_paths)
				if haskey(rhs_value.path_measurements, path_idx)
					qubits = rhs_value.path_measurements[path_idx]
					for qubit in qubits
						# Direct lookup in the measurements dictionary
						if haskey(sim.measurements[path_idx], qubit)
							outcome = sim.measurements[path_idx][qubit]
							set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, Any, outcome, false))
						end
					end
				end
			end
		else
			# Standard variable assignment
			for path_idx in sim.active_paths
				set_variable!(sim, path_idx, var_name, rhs_value)
			end
		end
	# This is for defining a specific bit in a register
	elseif head(lhs) == :indexed_identifier
		# Handle indexed assignment
		inds = _evolve_branched_ast(sim, lhs.args[2]) .+ 1  # Convert to 1-based indexing

		if is_measurement
			# For each active path, update the variable at the specified index
			for (path_idx_index, path_idx) in enumerate(sim.active_paths)
				if haskey(rhs_value.path_measurements, path_idx)
					qubits = rhs_value.path_measurements[path_idx]
					for qubit in qubits
						# Direct lookup in the measurements dictionary
						if haskey(sim.measurements[path_idx], qubit)
							outcome = sim.measurements[path_idx][qubit]
							
							# Get current value
							current_val = get_variable(sim, path_idx, var_name, nothing)

							# If the variable doesn't exist or isn't a vector, initialize it
							if isnothing(current_val) || !(current_val isa Vector)
								# Create a new vector with enough elements
								max_idx = maximum(inds)
								current_val = zeros(Int, max_idx)
							end

							# Update at index
							if current_val isa Vector
								current_val[inds] = outcome
								set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, Any, current_val, false))
							end
						end
					end
				end
			end
		else
			# Standard indexed assignment
			# For each active path, update the variable
			for path_idx in sim.active_paths
				# Get current value
				current_val = get_variable(sim, path_idx, var_name, nothing)

				# If the variable doesn't exist or isn't a vector, initialize it
				if isnothing(current_val) || !(current_val isa Vector)
					# Create a new vector with enough elements
					max_idx = maximum(inds)
					current_val = zeros(Int, max_idx)
				end

				# Update at index
				if current_val isa Vector
					current_val[inds] = rhs_value
					set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, Any, current_val, false))
				end
			end
		end
	end
end

"""
	_handle_classical_declaration(sim::BranchedSimulator, expr::QasmExpression)

Handle a classical variable declaration in the branched simulation model.
"""
function _handle_classical_declaration(sim::BranchedSimulator, expr::QasmExpression)
	if head(expr.args[2]) == :identifier
		var_name = Quasar.name(expr.args[2])
		var_type = head(expr.args[1]) == :classical_type ? expr.args[1].args[1] : expr.args[1]

		# Initialize variable for each active path
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, var_type, nothing, false))
		end
	elseif head(expr.args[2]) == :classical_assignment
		var_name = Quasar.name(expr.args[2].args[1].args[2])

		# Initialize variable for each active path and then handle the assignment
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, var_name, undef)
		end

		_handle_classical_assignment(sim, expr.args[2])
	end
end

"""
	_handle_const_declaration(sim::BranchedSimulator, expr::QasmExpression)

Handle a constant declaration in the branched simulation model.
"""
function _handle_const_declaration(sim::BranchedSimulator, expr::QasmExpression)
	if head(expr.args[2]) == :classical_assignment
		var_name = Quasar.name(expr.args[2].args[1].args[2])
		var_type = expr.args[1].args[1]

		# Initialize variable for each active path
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, var_name, undef)
		end

		# Handle the assignment
		_handle_classical_assignment(sim, expr.args[2])

		# Mark as constant in each path
		for path_idx in sim.active_paths
			if haskey(sim.variables[path_idx], var_name)
				sim.variables[path_idx][var_name].is_const = true
			end
		end
	else
		throw(QasmVisitorError("const declaration must assign an initial value."))
	end
end

"""
	_handle_qubit_declaration(sim::BranchedSimulator, expr::QasmExpression)

Handle a qubit declaration in the branched simulation model.
Dynamically resizes the simulator if needed to accommodate new qubits.
Uses the new Vector{Qubit} structure for qubit registers.
"""
function _handle_qubit_declaration(sim::BranchedSimulator, expr::QasmExpression)
	qubit_name = Quasar.name(expr)
	qubit_size = _evolve_branched_ast(sim, expr.args[2])

	for path_idx in sim.active_paths
		# Calculate the current qubit count
		current_qubit_count = sim.n_qubits[path_idx]

		# Create the qubit objects for this register
		qubit_objects = Vector{Qubit}()
		for i in 0:(qubit_size-1)
			qubit_idx = current_qubit_count + i
			qubit_obj = Qubit(qubit_idx, "$qubit_name[$i]", false)
			push!(qubit_objects, qubit_obj)
		end

		# Ensure qubit_mapping exists for this path
		if path_idx > length(sim.qubit_mapping)
			throw("Every simulation path must have a qubit mapping present")
		end

		# Store the entire register as a Vector{Qubit}
		sim.qubit_mapping[path_idx][qubit_name] = qubit_objects

		# For backward compatibility, also map individual qubits
		for (i, qubit_obj) in enumerate(qubit_objects)
			sim.qubit_mapping[path_idx]["$qubit_name[$(i-1)]"] = qubit_obj
		end

		new_qubit_count = current_qubit_count + qubit_size

		# Increment state size for each active path accordingly
		resize_simulator!(sim, path_idx, new_qubit_count)
	end
end

"""
	_handle_gate_modifiers(sim::BranchedSimulator, expr::QasmExpression)

Handle gate modifiers in the branched simulation model.
"""
function _handle_gate_modifiers(sim::BranchedSimulator, expr::QasmExpression)
	# Process gate modifiers
	mods = QasmExpression(:modifiers)
	mod_expr, inner = evaluate_modifiers(sim, expr)
	push!(mods, mod_expr)

	while head(inner) != :gate_call # done
		mod_expr, inner = evaluate_modifiers(sim, inner)
		push!(mods, mod_expr)
	end

	push!(inner, mods)
	_evolve_branched_ast(sim, inner)
end

"""
	_handle_gate_call(sim::BranchedSimulator, expr::QasmExpression)

Handle a gate call in the branched simulation model.
Instead of applying gates directly, store them in the instruction sequence.
When a measurement is needed, all stored instructions will be applied.
"""
function _handle_gate_call(sim::BranchedSimulator, expr::QasmExpression)
	gate_name = Quasar.name(expr)

	# Extract gate parameters
	params = []
	if length(expr.args) > 1 && !isnothing(expr.args[2])
		params = [_evolve_branched_ast(sim, param)
				  for param in convert(Vector{QasmExpression}, expr.args[2].args)]
	end

	# Extract target qubits for each path
	targets_by_path = []
	if length(expr.args) > 2 && !isnothing(expr.args[3])
		targets_by_path = evaluate_qubits(sim, expr.args[3])
	end

	# Get gate operator
	gate_op = nothing

	# Check if it's a built-in gate
	if haskey(sim.gate_defs, gate_name)
		gate_op = sim.gate_defs[gate_name](params...)
	else
		error("Gate $gate_name not defined!")
	end

	# Store gate instruction for each active path
	for (path_idx_index, path_idx) in enumerate(sim.active_paths)
		# Get targets for this specific path
		targets = isempty(targets_by_path) ? [] : targets_by_path[path_idx_index]
		
		# Check if any target qubit has been measured and removed
		mapped_targets = Int[]

		for qubit in targets
			# Check if this qubit was previously measured
			if qubit.measured
				# Expand the state to reincorporate the qubit
				expand_state!(sim, path_idx, qubit.index)
			end

			# Map the qubit to its current position in the state
			push!(mapped_targets, get_qubit_index(sim, path_idx, qubit))
		end

		# Create an instruction and add it to the sequence
		instruction = Instruction(gate_op, mapped_targets)
		push!(sim.instruction_sequences[path_idx], instruction)
	end

	return nothing
end

"""
	_handle_gate_definition(sim::BranchedSimulator, expr::QasmExpression)

Handle a gate definition in the branched simulation model.
"""
function _handle_gate_definition(sim::BranchedSimulator, expr::QasmExpression)
	# Register gate definition directly
	gate_name = Quasar.name(expr)
	gate_args = expr.args[2]
	gate_targets = expr.args[3]
	gate_body = expr.args[4]

	# Store the gate definition in the simulator
	argument_names = !isempty(gate_args.args) ?
					 [arg.args[1] for arg in gate_args.args[1]] : String[]
	qubit_targets = [Quasar.name(q) for q in gate_targets.args[1]]

	sim.gate_defs[gate_name] = GateDefinition(gate_name, argument_names, qubit_targets, gate_body)
end

"""
	_handle_function_call(sim::BranchedSimulator, expr::QasmExpression)

Handle a function call in the branched simulation model.
"""
function _handle_function_call(sim::BranchedSimulator, expr::QasmExpression)
	function_name = Quasar.name(expr)

	# Check if it's a built-in function
	if haskey(sim.function_defs, function_name)
		concrete_arguments = []
		if length(expr.args) > 1 && !isnothing(expr.args[2])
			concrete_arguments = [_evolve_branched_ast(sim, arg)
								  for arg in convert(Vector{QasmExpression}, expr.args[2].args)]
		end

		# Call the function
		return_val = sim.function_defs[function_name](concrete_arguments...)
		return return_val
	else
		error("Function $function_name not defined!")
	end
end

"""
	_handle_function_definition(sim::BranchedSimulator, expr::QasmExpression)

Handle a function definition in the branched simulation model.
"""
function _handle_function_definition(sim::BranchedSimulator, expr::QasmExpression)
	# Register function definition directly
	function_name = expr.args[1].args[1]
	function_args = expr.args[2]
	function_return_type = expr.args[3]
	function_body = expr.args[4]

	sim.function_defs[function_name] = FunctionDefinition(
		function_name, function_args, function_body, function_return_type)
end
