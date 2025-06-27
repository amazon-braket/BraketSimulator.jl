using Quasar: QasmExpression, head, ClassicalVariable, GateDefinition, FunctionDefinition, SizedBitVector

# Define a type to represent measurement values for each simulation path
# path_outcomes takes path indices as keys with the measurement results
# as dictionaries where the key is the qubit name and value is the qubit
# measurement
struct MeasurementValues
	path_outcomes::Dict{Int, Dict{String, Int}}
end

###############################################
# Helper functions for evaluating expressions #
###############################################

"""
	evaluate_qubits(sim::BranchedSimulatorOperators, qubit_expr::QasmExpression) -> Vector{Int}

Evaluate qubit expressions to get qubit indices. This replaces the visitor's evaluate_qubits function.
Returns a list of qubit indices in the state.
The same indices are used for all active simulations since we can't define qubits in a local scope.

As of now, this only works with registers and classical bits. Functionality for any classical variable
still needs to be implemented.
"""
function evaluate_qubits(sim::BranchedSimulatorOperators, qubit_expr::QasmExpression)
	# TODO: Update function to include functionality for any classical variable
	all_qubits = Int[]

	# Only process for active paths
	if isempty(sim.active_paths)
		return all_qubits
	end

	# Process the expression to get qubit indices
	expr_type = head(qubit_expr)

	if expr_type == :identifier
		qubit_name = Quasar.name(qubit_expr)

		# Check if the qubit name exists directly in the mapping
		if haskey(sim.qubit_mapping, qubit_name)
			qubit_idx = sim.qubit_mapping[qubit_name]
			push!(all_qubits, qubit_idx)
		else
			# Check if it's a register by looking for indexed qubits
			register_qubits = []
			for (key, value) in sim.qubit_mapping
				if startswith(key, "$qubit_name[")
					push!(register_qubits, value)
				end
			end

			if !isempty(register_qubits)
				append!(all_qubits, register_qubits)
			else
				error("Missing qubit '$qubit_name'")
			end
		end
	elseif expr_type == :indexed_identifier
		qubit_name = Quasar.name(qubit_expr)

		# Get the index
		qubit_ix = _evolve_branched_ast_operators(sim, qubit_expr.args[2])

		# Handle array indexing
		if qubit_ix isa Vector
			for ix in qubit_ix
				indexed_name = "$qubit_name[$ix]"
				haskey(sim.qubit_mapping, indexed_name) || error("Missing qubit '$indexed_name'")
				qubit_idx = sim.qubit_mapping[indexed_name]
				push!(all_qubits, qubit_idx)
			end
		else
			indexed_name = "$qubit_name[$qubit_ix]"
			haskey(sim.qubit_mapping, indexed_name) || error("Missing qubit '$indexed_name'")
			qubit_idx = sim.qubit_mapping[indexed_name]
			push!(all_qubits, qubit_idx)
		end
	elseif expr_type == :array_literal
		# For array literals, recursively evaluate each element
		for element in qubit_expr.args
			append!(all_qubits, evaluate_qubits(sim, element))
		end
	elseif expr_type == :hw_qubit
		# Hardware qubit reference (e.g., $0, $1)
		qubit_idx = parse(Int, replace(qubit_expr.args[1], "\$" => ""))
		push!(all_qubits, qubit_idx)
	end

	return all_qubits
end

"""
	evaluate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)

Evaluates gate modifiers.
"""
function evaluate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Review and edit this function
	if head(expr) == :power_mod
		pow_expr = QasmExpression(:pow, _evolve_branched_ast_operators(sim, expr.args[1]))
		return (pow_expr, expr.args[2])
	elseif head(expr) == :inverse_mod
		return (QasmExpression(:inv), expr.args[1])
	elseif head(expr) ∈ (:control_mod, :negctrl_mod)
		has_argument = length(expr.args) > 1
		if has_argument
			arg_val = _evolve_branched_ast_operators(sim, first(expr.args))
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
	evolve_branched_operators(simulator::AbstractSimulator, program::QasmExpression, inputs::Dict{String, <:Any}) -> BranchedSimulatorOperators

Evolve a quantum program using a branched approach that handles measurements and control flow.

Takes in an AbstractSimulator object and generates a branched simulator object to wrap around it in order
to perform the MCM. Allows it to be generalizable if you have any simulator as long as it is a subtype of AbstractSimulator.
"""
function evolve_branched_operators(simulator::AbstractSimulator, program::QasmExpression, inputs::Dict{String, <:Any})
	# Create a branched simulator with integrated visitor functionality
	branched_sim = BranchedSimulatorOperators(simulator; inputs = inputs)

	# Process the AST
	_evolve_branched_ast_operators(branched_sim, program)

	return branched_sim
end


_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::Number) = i
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::String) = i
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::BitVector) = i
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::NTuple{N, <:Number}) where {N} = i
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::Vector{<:Number}) = i

"""
	_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, expr::QasmExpression)

Process an AST node in the branched simulation model. AST taken from Quasar.jl
parse function.
"""
function _evolve_branched_ast_operators(sim::BranchedSimulatorOperators, expr::QasmExpression)
	expr_type = head(expr)

	######################################
	# Organizational/Unimplemented nodes #
	######################################

	# Program node - process each child node in sequence
	if expr_type == :program
		for child_expr in expr.args
			head(child_expr) == :end && return
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && return  # No active paths left
		end

		# Scope node - process each child node in sequence
	elseif expr_type == :scope
		for child_expr in expr.args
			child_head = head(child_expr)
			if child_head in (:end, :continue, :break)
				return
			end
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && return  # No active paths left
		end

	elseif expr_type == :version
		# Version node - ignore
		return

		# TODO: Have not looked at this code, need to check it and see if it is implemented
	elseif expr_type == :reset
		# Reset operation - reset qubit to |0⟩ state
		qubit_indices = evaluate_qubits(sim, expr.args[1])
		for path_idx in sim.active_paths
			for qubit_idx in qubit_indices
				# Apply reset to the qubit's index in the state
				_apply_reset!(path_idx, qubit_idx)
			end
		end

	elseif expr_type == :barrier || expr_type == :delay || expr_type == :pragma
		# Barrier/delay/pragma operation - no effect on simulation
		return

	elseif expr_type in [:stretch, :duration, :box, :output]
		# Unsupported operations
		error("Simulator doesn't support $expr_type operations")

		################################################
		# Variable declarations and control flow nodes #
		################################################
	elseif expr_type == :input
		# Input parameter declaration for parameterized functions
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

	elseif expr_type == :for
		_handle_for_loop(sim, expr)

	elseif expr_type == :switch
		_handle_switch_statement(sim, expr)

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
		var = get_variable(sim, path_idx, id_name)

		# This is because we store the value of singular bits as a bit register with one bit
		if var.type == Quasar.SizedBitVector && length(var.val) == 1
			return var.val[1]
		end

		if !isnothing(var.val)
			return var.val
		end

		# Then check qubit_mapping
		if haskey(sim.qubit_mapping, id_name)
			return sim.qubit_mapping[id_name]
		end

		# If we get here, the identifier is not defined
		error("No identifier $id_name defined")

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
		index = _evolve_branched_ast_operators(sim, expr.args[2])

		# Get the variable
		var = get_variable(sim, path_idx, identifier_name)

		if !isnothing(var)
			# Check if it's a bit array (SizedBitVector)
			if var.type == Quasar.SizedBitVector
				# Convert to 1-based indexing for Julia arrays
				julia_index = index + 1

				if julia_index <= length(var.val)
					# Access the bit directly from the boolean array
					return var.val[julia_index]
				end
			elseif var.type isa Vector
				# For other vector types
				# Convert to 1-based indexing for Julia arrays
				flat_ix = (index isa Vector) ? index : [index]
				flat_ix = flat_ix .+ 1

				return var.val[only(flat_ix)]
			end
		end

		# Then check qubit_mapping
		return evaluate_qubits(sim, expr)

	elseif expr_type == :if
		_handle_conditional(sim, expr)

	elseif expr_type == :while
		_handle_while_loop(sim, expr)

	elseif expr_type == :classical_assignment
		_handle_classical_assignment(sim, expr)

	elseif expr_type == :classical_declaration
		_handle_classical_declaration(sim, expr)

	elseif expr_type == :const_declaration
		_handle_const_declaration(sim, expr)

	elseif expr_type == :qubit_declaration
		_handle_qubit_declaration(sim, expr)

	elseif expr_type ∈ [:power_mod, :inverse_mod, :control_mod, :negctrl_mod]
		_handle_gate_modifiers(sim, expr)

	elseif expr_type == :gate_call
		_handle_gate_call(sim, expr)

	elseif expr_type == :gate_definition
		_handle_gate_definition(sim, expr)

	elseif expr_type == :measure
		return _handle_measurement(sim, expr)

	elseif expr_type == :function_call
		_handle_function_call(sim, expr)

	elseif expr_type == :function_definition
		_handle_function_definition(sim, expr)

	elseif expr_type ∈ [:integer_literal, :float_literal, :string_literal, :complex_literal, :irrational_literal, :boolean_literal, :duration_literal]
		# Return the literal value
		return expr.args[1]

	elseif expr_type == :array_literal
		# Evaluate array elements
		return [_evolve_branched_ast_operators(sim, element) for element in expr.args]

	elseif expr_type == :range
		# Evaluate range parameters
		start = _evolve_branched_ast_operators(sim, expr.args[1])
		stop = _evolve_branched_ast_operators(sim, expr.args[2])
		step = length(expr.args) > 2 ? _evolve_branched_ast_operators(sim, expr.args[3]) : 1
		return collect(start:step:stop)

	elseif expr_type == :binary_op
		# Binary operation
		op = expr.args[1]
		lhs = _evolve_branched_ast_operators(sim, expr.args[2])
		rhs = _evolve_branched_ast_operators(sim, expr.args[3])

		# Evaluate based on operator type
		return evaluate_binary_op(op, lhs, rhs)

	elseif expr_type == :unary_op
		# Unary operation
		op = expr.args[1]
		arg = _evolve_branched_ast_operators(sim, expr.args[2])

		# Evaluate using existing function
		return evaluate_unary_op(op, arg)

	elseif expr_type == :cast
		casting_to = expr.args[1].args[1]
		value = _evolve_branched_ast_operators(sim, expr.args[2])

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
	_handle_measurement(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a measurement operation, creating branches for different outcomes.
Keeps the full state size and sets probabilities to 0 for states that don't match the measurement outcome.
Returns a MeasurementValues object containing the measurement outcomes for each path.

We assume that the qubit indices to measure is identical among all of the paths because
of openqasm restrictions regarding scope.
"""
function _handle_measurement(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# Get qubit indices to measure
	qubit_indices = evaluate_qubits(sim, expr.args[1])

	new_active_paths = Int[]

	# Structure to store measurement outcomes for each path
	path_outcomes = Dict{Int, Dict{String, Int}}()

	for path_idx in sim.active_paths
		# Everytime a measurement occurs, the resulting path indices are added to the following list
		paths_to_measure = [path_idx]

		path_outcomes[path_idx] = Dict{String, Int}()

		# For every active path, measure each qubit in the qubit indices
		for qubit_idx in qubit_indices
			qubit_name = ""
			for (name, idx) in sim.qubit_mapping
				if idx == qubit_idx
					qubit_name = name
					break
				end
			end

			# This inner for loop iterates through all sub paths created by a given path undergoing a measurement
			current_paths = copy(paths_to_measure)
			for idx in current_paths
				# Check if measurement caused a new path (both measurement outcome posible)
				path_split = measure_qubit(sim, idx, qubit_idx, qubit_name)

				outcome = sim.measurements[idx][qubit_name][end]
				path_outcomes[idx][qubit_name] = outcome

				# Create and store a measure operator with the result
				measure_op = Measure(outcome)
				instruction = Instruction(measure_op, [qubit_idx])
				push!(sim.instruction_sequences[idx], instruction)

				if path_split
					# A new path was created
					added_idx = length(sim.instruction_sequences) # Since we add the index at the end of the states vector
					push!(paths_to_measure, added_idx)

					# Copy all previous outcomes to the new path
					path_outcomes[added_idx] = copy(path_outcomes[idx])

					# Get and update the measurement outcome for the new path and qubit
					new_outcome = 1 - outcome
					path_outcomes[added_idx][qubit_name] = new_outcome

					# Create and store a measure operator with the result for the new path
					new_measure_op = Measure(new_outcome)
					new_instruction = Instruction(new_measure_op, [qubit_idx])
					push!(sim.instruction_sequences[added_idx], new_instruction)
				end
			end
		end

		# Add all the new subpaths to the active paths list
		append!(new_active_paths, paths_to_measure)
	end

	# Update active paths all at once
	sim.active_paths = new_active_paths

	return MeasurementValues(path_outcomes)
end

"""
	_handle_conditional(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a conditional statement, filtering paths based on the condition.
"""
function _handle_conditional(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# Find if there's an else branch
	has_else = findfirst(e->head(e) == :else, convert(Vector{QasmExpression}, expr.args))
	last_expr = !isnothing(has_else) ? length(expr.args) - 1 : length(expr.args)

	# Evaluate condition for each active path
	true_paths = Int[]
	false_paths = Int[]

	for path_idx in sim.active_paths
		condition_value = _evaluate_condition(sim, path_idx, expr.args[1])

		if condition_value
			push!(true_paths, path_idx)
		else
			push!(false_paths, path_idx)
		end
	end

	# Process if branch for true paths
	if !isempty(true_paths)
		original_active_paths = copy(sim.active_paths)
		sim.active_paths = true_paths

		# Process if branch
		for child_expr in expr.args[2:last_expr]
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && break
		end

		# Restore active paths that were not processed in if branch
		append!(sim.active_paths, setdiff(original_active_paths, true_paths))
	end

	# Process else branch for false paths
	if !isnothing(has_else) && !isempty(false_paths)
		original_active_paths = copy(sim.active_paths)
		sim.active_paths = false_paths

		# Process else branch
		for child_expr in expr.args[has_else].args
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && break
		end

		# Restore active paths that were not processed in else branch
		append!(sim.active_paths, setdiff(original_active_paths, false_paths))
	end
end

"""
	_evaluate_condition(sim::BranchedSimulatorOperators, path_idx::Int, condition_expr::QasmExpression)

Evaluate a condition expression in the context of a specific path.
This handles path-specific variables and measurement results.
"""
function _evaluate_condition(sim::BranchedSimulatorOperators, path_idx::Int, condition_expr::QasmExpression)
	# Save the current active paths
	original_active_paths = copy(sim.active_paths)

	sim.active_paths = [path_idx]
	result = _evolve_branched_ast_operators(sim, condition_expr)

	# Restore active paths
	sim.active_paths = original_active_paths

	return result isa Bool ? result : result > 0
end

"""
	_handle_for_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a for loop in the branched simulation model.
"""
function _handle_for_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check for loop implementation
	for_loop = convert(Vector{QasmExpression}, expr.args)
	loop_variable_type = for_loop[1].args[1]
	loop_variable_name = for_loop[2].args[1]::String

	# Evaluate loop range
	loop_variable_values = _evolve_branched_ast_operators(sim, for_loop[3])
	loop_body = for_loop[4]::QasmExpression

	# For each value in the range
	for loop_value in loop_variable_values
		# Set the variable for each active path
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, loop_variable_name, ClassicalVariable(loop_variable_name, loop_variable_type, loop_value, false))
		end

		# Process loop body
		_evolve_branched_ast_operators(sim, loop_body)

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
	_handle_switch_statement(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a switch statement in the branched simulation model. For more information, visit
https://openqasm.com/language/classical.html#the-switch-statement
"""
function _handle_switch_statement(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check switch statement implementation
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
		case_val = _evolve_branched_ast_operators(sim, expr.args[1])

		# Find matching case for this path
		path_matched = false

		for (case_idx, case) in enumerate(all_cases)
			if head(case) == :case
				# Evaluate case value for this path
				sim.active_paths = [path_idx]
				case_value_result = _evolve_branched_ast_operators(sim, case.args[1])

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
			_evolve_branched_ast_operators(sim, child_expr)
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
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && break  # No active paths left for default
		end

		# Add surviving paths to the collection
		append!(surviving_paths, sim.active_paths)
	elseif !isempty(default_paths) && isnothing(default_idx)
		# Just do a no operation, and add the default paths onto the active paths
		append!(surviving_paths, default_paths)
	end

	# Update active paths to include all paths that survived their respective cases
	sim.active_paths = surviving_paths

	# If no paths survived, restore original paths
	if isempty(sim.active_paths)
		sim.active_paths = all_active_paths
	end
end

"""
	_handle_while_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a while loop in the branched simulation model.
"""
function _handle_while_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check while loop implementation
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
			condition_value = _evolve_branched_ast_operators(sim, expr.args[1]) > 0

			if condition_value
				push!(new_continue_paths, path_idx)
			end
		end

		# If no paths should continue, break
		isempty(new_continue_paths) && break

		# Set active paths to only those that should continue
		sim.active_paths = new_continue_paths

		# Process loop body
		_evolve_branched_ast_operators(sim, loop_body)

		# Update continue_paths for next iteration
		continue_paths = copy(sim.active_paths)

		# Restore paths that didn't enter the loop
		append!(sim.active_paths, setdiff(original_active_paths, new_continue_paths))
	end
end

"""
	_handle_classical_assignment(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a classical variable assignment in the branched simulation model.
"""
function _handle_classical_assignment(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Update this to include all classical variable types
	op = expr.args[1].args[1]
	lhs = expr.args[1].args[2]
	rhs = expr.args[1].args[3]

	var_name = Quasar.name(lhs)

	new_active_paths = []

	# Iterate through each active path
	for path_idx in copy(sim.active_paths)
		# Set this path as the only active one temporarily for evaluation
		sim.active_paths = [path_idx]

		# Evaluate the right-hand side for this specific path
		# This will handle any measurements or complex expressions
		rhs_value = _evolve_branched_ast_operators(sim, rhs)

		# Check if the right-hand side is a measurement result
		is_measurement = rhs_value isa MeasurementValues

		# Get all paths that were created during RHS evaluation
		current_paths = copy(sim.active_paths)

		# Process each path (original and any new ones created during RHS evaluation)
		for current_path_idx in current_paths
			# This is for assigning to a bit register or any classical variable
			if head(lhs) == :identifier
				if is_measurement
					# Process measurement outcomes for this path
					if haskey(rhs_value.path_outcomes, current_path_idx)
						outcomes = rhs_value.path_outcomes[current_path_idx]

						# Get the existing variable
						var = get(sim.variables[current_path_idx], var_name, nothing)

						if !isnothing(var) && !isempty(outcomes)
							# Sort outcomes by qubit name for consistent ordering
							sorted_outcomes = sort(collect(outcomes))

							# Extract bit values and reverse them to use little endian notation
							# This matches the array access pattern which is little endian
							bit_values = reverse([outcome for (_, outcome) in sorted_outcomes])

							# The register size should already be set during initialization
							# We're just updating the boolean array with the measurement values
							bit_array = var.val
							if bit_array isa Vector{Int} && length(bit_array) == length(bit_values)
								# Update the boolean array directly
								var.val = bit_values
							else
								error("Assignment must be to a register of the right size")
							end
						end
					end
				else
					# Standard variable assignment
					# Get the existing variable
					var = get(sim.variables[current_path_idx], var_name, nothing)

					if !isnothing(var)
						# Update the existing variable's value, preserving its type
						var.val = rhs_value
					end
				end
			# This is for assigning to a specific bit in a register
			elseif head(lhs) == :indexed_identifier
				# Handle indexed assignment
				index = _evolve_branched_ast_operators(sim, lhs.args[2])  # Get the 0-based index

				if is_measurement
					# Process measurement for this path
					if haskey(rhs_value.path_outcomes, current_path_idx)
						outcomes = rhs_value.path_outcomes[current_path_idx]

						# Get the existing register variable
						register_var = get(sim.variables[current_path_idx], var_name, nothing)

						if !isnothing(register_var) && !isempty(outcomes)
							# Use the first measurement outcome for the specified index
							_, outcome = first(outcomes)

							# Update the bit directly in the boolean array
							bit_array = register_var.val
							if bit_array isa Vector{Int} && index + 1 <= length(bit_array)
								# Convert to 1-based indexing for Julia arrays
								julia_index = index + 1
								bit_array[julia_index] = outcome
							end
						end
					end
				else
					# Standard indexed assignment
					# Update the bit directly in the boolean array
					# Get the register variable
					register_var = get(sim.variables[current_path_idx], var_name, nothing)

					if !isnothing(register_var)
						bit_array = register_var.val
						if bit_array isa Vector{Int} && index + 1 <= length(bit_array)
							# Convert to 1-based indexing for Julia arrays
							julia_index = index + 1
							bit_array[julia_index] = rhs_value
						end
					end
				end
			end
		end
		# Add all the active paths generated in the new_active_paths variable
		append!(new_active_paths, sim.active_paths)
	end
	sim.active_paths = copy(new_active_paths)
end

"""
	_handle_classical_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a classical variable declaration in the branched simulation model.
"""
function _handle_classical_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Update this to include all classical variable types
	# This is for a register/bit
	if head(expr.args[2]) == :identifier
		var_name = Quasar.name(expr.args[2])
		var_type = head(expr.args[1]) == :classical_type ? expr.args[1].args[1] : expr.args[1]

		# Check if this is a bit register declaration
		if typeof(var_type) == Quasar.SizedBitVector
			# Get the size of the register
			register_size = _evolve_branched_ast_operators(sim, var_type.size)

			# If the size is -1, then it is a single bit, not a register
			if register_size == -1
				# Initialize a single bit for each active path
				for path_idx in sim.active_paths
					# Initialize with a single boolean value
					set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, typeof(var_type), [1], false))
				end
			else
				# Initialize the register for each active path
				for path_idx in sim.active_paths
					# Initialize the register with a boolean array
					set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, typeof(var_type), [0 for _ in 1:register_size], false))
				end
			end
		else
			# Initialize variable for each active path
			for path_idx in sim.active_paths
				set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, var_type, nothing, false))
			end
		end
		# This for all other classical variables
	elseif head(expr.args[2]) == :classical_assignment
		var_name = Quasar.name(expr.args[2].args[1].args[2])
		var_type = expr.args[1].args[1]

		# Initialize variable for each active path and then handle the assignment
		for path_idx in sim.active_paths
			set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, var_type, nothing, false))
		end

		_handle_classical_assignment(sim, expr.args[2])
	end
end

"""
	_handle_const_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a constant declaration in the branched simulation model.
"""
function _handle_const_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Update this to include all classical variable types
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
	_handle_qubit_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a qubit declaration in the branched simulation model.
Maps qubit names to their indices in the state vector.
"""
function _handle_qubit_declaration(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Update this to include qubit declaration with -1 index (individual qubit)
	qubit_name = Quasar.name(expr)
	qubit_size = _evolve_branched_ast_operators(sim, expr.args[2])
	current_qubit_count = sim.n_qubits

	# Map qubit indices
	for i in 0:(qubit_size-1)
		qubit_idx = current_qubit_count + i
		indexed_name = "$qubit_name[$i]"
		sim.qubit_mapping[indexed_name] = qubit_idx
	end

	# Store the base name as well if it's a single qubit
	if qubit_size == -1
		sim.qubit_mapping[qubit_name] = current_qubit_count
	end

	sim.n_qubits = current_qubit_count + qubit_size
end

"""
	_handle_gate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle gate modifiers in the branched simulation model.
"""
function _handle_gate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check this function
	# Process gate modifiers
	mods = QasmExpression(:modifiers)
	mod_expr, inner = evaluate_modifiers(sim, expr)
	push!(mods, mod_expr)

	while head(inner) != :gate_call # done
		mod_expr, inner = evaluate_modifiers(sim, inner)
		push!(mods, mod_expr)
	end

	push!(inner, mods)
	_evolve_branched_ast_operators(sim, inner)
end

"""
	_handle_gate_call(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a gate call in the branched simulation model.
Instead of applying gates directly, store them in the instruction sequence.
When a measurement is needed, all stored instructions will be applied.
"""
function _handle_gate_call(sim::BranchedSimulatorOperators, expr::QasmExpression)
	gate_name = Quasar.name(expr)

	# Extract gate parameters
	params = []
	if length(expr.args) > 1 && !isnothing(expr.args[2])
		params = [_evolve_branched_ast_operators(sim, param)
				  for param in convert(Vector{QasmExpression}, expr.args[2].args)]
	end

	# Extract target qubit indices
	target_indices = []
	if length(expr.args) > 2 && !isnothing(expr.args[3])
		if head(expr.args[3]) == :qubit_targets
			for target in expr.args[3].args
				qubit_indices = evaluate_qubits(sim, target)
				append!(target_indices, qubit_indices)
			end
		else
			# Handle the case where it's not a :qubit_targets expression
			target_indices = evaluate_qubits(sim, expr.args[3])
		end
	end

	# Get gate operator
	gate_op = nothing
	if haskey(sim.gate_defs, gate_name)
		gate_def = sim.gate_defs[gate_name]
		instruction = gate_def.body
		gate_op = StructTypes.constructfrom(QuantumOperator, instruction)
	else
		error("Gate $gate_name not defined!")
	end

	# Store gate instruction for each active path
	for path_idx in sim.active_paths
		instruction = Instruction(gate_op, target_indices)
		push!(sim.instruction_sequences[path_idx], instruction)
	end

	return nothing
end

"""
	_handle_gate_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a gate definition in the branched simulation model.
"""
function _handle_gate_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check this function implementation
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
	_handle_function_call(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a function call in the branched simulation model.
"""
function _handle_function_call(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check this function implementation
	function_name = Quasar.name(expr)

	# Check if it's a built-in function
	if haskey(sim.function_defs, function_name)
		concrete_arguments = []
		if length(expr.args) > 1 && !isnothing(expr.args[2])
			concrete_arguments = [_evolve_branched_ast_operators(sim, arg)
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
	_handle_function_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a function definition in the branched simulation model.
"""
function _handle_function_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# TODO: Check this function implementation
	# Register function definition directly
	function_name = expr.args[1].args[1]
	function_args = expr.args[2]
	function_return_type = expr.args[3]
	function_body = expr.args[4]

	sim.function_defs[function_name] = FunctionDefinition(
		function_name, function_args, function_body, function_return_type)
end
