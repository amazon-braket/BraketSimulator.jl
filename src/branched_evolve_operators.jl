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
	evaluate_qubits(sim::BranchedSimulatorOperators, qubit_expr::QasmExpression) -> Dict{Int, Vector{Int}}

Evaluate qubit expressions to get qubit indices. This replaces the visitor's evaluate_qubits function.
Returns a dictionary mapping path indices to lists of qubit indices in the state.
The same indices are used for all active simulations since we can't define qubits in a local scope.

As of now, this only works with registers and classical bits. Functionality for any classical variable
still needs to be implemented.
"""
function evaluate_qubits(sim::BranchedSimulatorOperators, qubit_expr::QasmExpression)
	# TODO: Update function to include functionality for any classical variable
	# Create a dictionary to store results for each path
	results = Dict{Int, Vector{Int}}()

	# Only process for active paths
	if isempty(sim.active_paths)
		return results
	end

	# Process the expression to get qubit indices
	expr_type = head(qubit_expr)

	if expr_type == :identifier
		qubit_name = Quasar.name(qubit_expr)

		for path_idx in sim.active_paths
			all_qubits = Int[]

			# First check if it's a function parameter that refers to a qubit
			var = get_variable(sim, path_idx, qubit_name)

			if !isnothing(var) && var.type == :qubit_declaration && !isnothing(var.val)
				# If it's a qubit parameter, use the stored qubit name
				stored_qubit_idx = var.val
				if stored_qubit_idx isa Vector
					append!(all_qubits, stored_qubit_idx)
				else
					push!(all_qubits, stored_qubit_idx)
				end
				# Check if it's a qubit alias
			elseif !isnothing(var) && var.type == :qubit_alias && !isnothing(var.val)
				# If it's a qubit alias, use the stored qubit indices
				alias_qubits = var.val
				if alias_qubits isa Vector
					append!(all_qubits, alias_qubits)
				else
					push!(all_qubits, alias_qubits)
				end
				# Check if the qubit name exists directly in the mapping
			elseif haskey(sim.qubit_mapping, qubit_name)
				qubit_idx = sim.qubit_mapping[qubit_name]
				if qubit_idx isa Vector
					append!(all_qubits, qubit_idx)
				else
					push!(all_qubits, qubit_idx)
				end
			else
				# Check if it's a register by looking for indexed qubits
				register_qubits = []
				for (key, value) in sim.qubit_mapping
					if startswith(key, "$qubit_name[")
						push!(register_qubits, value)
					end
				end

				if !isempty(register_qubits)
					for qubit_idx in register_qubits
						push!(all_qubits, qubit_idx)
					end
				else
					error("Missing qubit '$qubit_name'")
				end
			end

			results[path_idx] = all_qubits
		end
	elseif expr_type == :indexed_identifier
		qubit_name = Quasar.name(qubit_expr)

		# Get the index
		qubit_ix = _evolve_branched_ast_operators(sim, qubit_expr.args[2])

		for path_idx in sim.active_paths
			all_qubits = Int[]

			# Get the index for this specific path
			path_qubit_ix = if isa(qubit_ix, Dict)
				# If it's a dictionary, get the value for this path
				if haskey(qubit_ix, path_idx)
					qubit_ix[path_idx]
				else
					# If this path isn't in the dictionary, use the first value
					first(values(qubit_ix))
				end
			else
				qubit_ix
			end

			# First check if it's a qubit alias in variables
			var = get_variable(sim, path_idx, qubit_name)
			if !isnothing(var) && var.type == :qubit_alias && !isnothing(var.val)
				# If the base name is an alias for a register, use the indexed qubit from the register
				register_qubits = var.val

				# Handle array indexing
				if path_qubit_ix isa Vector
					for ix in path_qubit_ix
						# Make sure the index is valid
						if ix >= 0 && ix < length(register_qubits)
							push!(all_qubits, register_qubits[ix+1])  # Convert to 1-based indexing
						else
							error("Index $ix out of bounds for register $qubit_name")
						end
					end
				else
					# Make sure the index is valid
					if path_qubit_ix >= 0 && path_qubit_ix < length(register_qubits)
						push!(all_qubits, register_qubits[path_qubit_ix+1])  # Convert to 1-based indexing
					else
						error("Index $path_qubit_ix out of bounds for register $qubit_name")
					end
				end
				# Then check if the base name is an alias for a register in qubit_mapping
			elseif haskey(sim.qubit_mapping, qubit_name) && sim.qubit_mapping[qubit_name] isa Vector
				# If the base name is an alias for a register, use the indexed qubit from the register
				register_qubits = sim.qubit_mapping[qubit_name]

				# Handle array indexing
				if path_qubit_ix isa Vector
					for ix in path_qubit_ix
						# Make sure the index is valid
						if ix >= 0 && ix < length(register_qubits)
							push!(all_qubits, register_qubits[ix+1])  # Convert to 1-based indexing
						else
							error("Index $ix out of bounds for register $qubit_name")
						end
					end
				else
					# Make sure the index is valid
					if path_qubit_ix >= 0 && path_qubit_ix < length(register_qubits)
						push!(all_qubits, register_qubits[path_qubit_ix+1])  # Convert to 1-based indexing
					else
						error("Index $path_qubit_ix out of bounds for register $qubit_name")
					end
				end
			else
				# Handle array indexing for regular qubits
				if path_qubit_ix isa Vector
					for ix in path_qubit_ix
						indexed_name = "$qubit_name[$ix]"
						haskey(sim.qubit_mapping, indexed_name) || error("Missing qubit '$indexed_name'")
						qubit_idx = sim.qubit_mapping[indexed_name]
						push!(all_qubits, qubit_idx)
					end
				else
					indexed_name = "$qubit_name[$path_qubit_ix]"
					haskey(sim.qubit_mapping, indexed_name) || error("Missing qubit '$indexed_name'")
					qubit_idx = sim.qubit_mapping[indexed_name]
					push!(all_qubits, qubit_idx)
				end
			end

			results[path_idx] = all_qubits
		end
	elseif expr_type == :array_literal
		# For array literals, recursively evaluate each element
		for path_idx in sim.active_paths
			all_qubits = Int[]

			for element in qubit_expr.args
				# Save the original active paths
				original_active_paths = copy(sim.active_paths)

				# Temporarily set this path as the only active one for evaluation
				sim.active_paths = [path_idx]

				# Get qubits for this specific path
				path_qubits = evaluate_qubits(sim, element)

				# Restore original active paths
				sim.active_paths = original_active_paths

				if !isempty(path_qubits) && haskey(path_qubits, path_idx)
					append!(all_qubits, path_qubits[path_idx])
				end
			end

			results[path_idx] = all_qubits
		end
	elseif expr_type == :hw_qubit
		# Hardware qubit reference (e.g., $0, $1)
		qubit_idx = parse(Int, replace(qubit_expr.args[1], "\$" => ""))

		for path_idx in sim.active_paths
			results[path_idx] = [qubit_idx]
		end
	end

	return results
end

"""
	evaluate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)

Evaluates gate modifiers. Always returns a dictionary mapping path indices to tuples of (modifier_expr, inner_expr).
"""
function evaluate_modifiers(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# Create a dictionary to store results for each path
	results = Dict{Int, Tuple{QasmExpression, Any}}()

	if head(expr) == :power_mod
		pow_val = _evolve_branched_ast_operators(sim, expr.args[1])

		# Check if pow_val is a dictionary (path-specific values)
		if isa(pow_val, Dict)
			# For each path, evaluate with path-specific value
			for (path_idx, path_val) in pow_val
				pow_expr = QasmExpression(:pow, path_val)
				results[path_idx] = (pow_expr, expr.args[2])
			end
		else
			# If not a dictionary, use the same value for all active paths
			pow_expr = QasmExpression(:pow, pow_val)
			for path_idx in sim.active_paths
				results[path_idx] = (pow_expr, expr.args[2])
			end
		end

		return results
	elseif head(expr) == :inverse_mod
		# Use the same inverse modifier for all active paths
		for path_idx in sim.active_paths
			results[path_idx] = (QasmExpression(:inv), expr.args[1])
		end

		return results
	elseif head(expr) ∈ (:control_mod, :negctrl_mod)
		has_argument = length(expr.args) > 1
		if has_argument
			arg_val = _evolve_branched_ast_operators(sim, first(expr.args))

			# Check if arg_val is a dictionary (path-specific values)
			if isa(arg_val, Dict)
				# For each path, evaluate with path-specific value
				for (path_idx, path_val) in arg_val
					isinteger(path_val) || error("Cannot apply non-integer ($path_val) number of controls or negcontrols.")
					true_inner = expr.args[2]
					inner = QasmExpression(head(expr), true_inner)

					local_arg_val = path_val
					while local_arg_val > 2
						inner = QasmExpression(head(expr), inner)
						local_arg_val -= 1
					end

					new_head = head(expr) == :control_mod ? :ctrl : :negctrl
					results[path_idx] = (QasmExpression(new_head), inner)
				end
			else
				# If not a dictionary, use the same value for all active paths
				isinteger(arg_val) || error("Cannot apply non-integer ($arg_val) number of controls or negcontrols.")
				true_inner = expr.args[2]
				inner = QasmExpression(head(expr), true_inner)

				while arg_val > 2
					inner = QasmExpression(head(expr), inner)
					arg_val -= 1
				end

				new_head = head(expr) == :control_mod ? :ctrl : :negctrl

				for path_idx in sim.active_paths
					results[path_idx] = (QasmExpression(new_head), inner)
				end
			end
		else
			inner = expr.args[1]
			new_head = head(expr) == :control_mod ? :ctrl : :negctrl

			for path_idx in sim.active_paths
				results[path_idx] = (QasmExpression(new_head), inner)
			end
		end

		return results
	else
		error("Unknown modifier type: $(head(expr))")
	end
end

"""
	_apply_modifiers(gate_op::QuantumOperator, modifiers::Vector, target_indices::Vector{Int}) -> QuantumOperator

Apply gate modifiers to a gate operator. This function takes a gate operator and a list of modifiers,
and returns a new gate operator with the modifiers applied.
"""
function _apply_modifiers(gate_op::QuantumOperator, modifiers::Vector, target_indices::Vector{Int})
	# Process each modifier in order
	for modifier in modifiers
		modifier_type = head(modifier)

		if modifier_type == :ctrl
			# Special handling for global phase gates
			if isa(gate_op, GPhase)
				if length(target_indices) < 1
					error("Control modifier requires at least 1 qubit for phase gates")
				end

				# Create a controlled version of the gate
				# For phase gates, the control is the only qubit
				bitvals = tuple(1)

				# Create the controlled gate
				gate_op = Control(gate_op, bitvals)
			else
				# For regular gates, we need at least 2 qubits (control + target)
				if length(target_indices) < 2
					error("Control modifier requires at least 2 qubits")
				end

				# Create a controlled version of the gate
				# The first qubit is the control, the rest are targets for the gate
				bitvals = tuple(1)

				# Create the controlled gate
				gate_op = Control(gate_op, bitvals)
			end

		elseif modifier_type == :negctrl
			# Special handling for global phase gates
			if isa(gate_op, PhaseShift)
				# For PhaseShift (gphase), we only need one control qubit
				if length(target_indices) < 1
					error("Negative control modifier requires at least 1 qubit for phase gates")
				end

				# Create a controlled version of the gate with control on |0⟩ state
				bitvals = tuple(0)

				# Create the negatively controlled gate
				gate_op = Control(gate_op, bitvals)
			else
				# Apply negative control modifier (control on |0⟩ state)
				if length(target_indices) < 2
					error("Negative control modifier requires at least 2 qubits")
				end

				# Create a controlled version of the gate with control on |0⟩ state
				bitvals = tuple(0)

				# Create the negatively controlled gate
				gate_op = Control(gate_op, bitvals)
			end

		elseif modifier_type == :pow
			# Apply power modifier
			if hasfield(typeof(gate_op), :pow_exponent)
				# If the gate has a pow_exponent field, update it
				power_value = modifier.args[1]
				gate_op.pow_exponent = power_value
			else
				# If the gate doesn't have a pow_exponent field, we can't apply the power modifier
				error("Power modifier not supported for gate type $(typeof(gate_op))")
			end

		elseif modifier_type == :inv
			# Apply inverse modifier
			if hasfield(typeof(gate_op), :pow_exponent)
				# If the gate has a pow_exponent field, negate it to invert
				gate_op.pow_exponent = -gate_op.pow_exponent
			else
				# If the gate doesn't have a pow_exponent field, we can't apply the inverse modifier
				error("Inverse modifier not supported for gate type $(typeof(gate_op))")
			end
		end
	end

	return gate_op
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


_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::Number) = Dict(path_idx => i for path_idx in sim.active_paths)
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::String) = Dict(path_idx => i for path_idx in sim.active_paths)
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::BitVector) = Dict(path_idx => i for path_idx in sim.active_paths)
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::NTuple{N, <:Number}) where {N} = Dict(path_idx => i for path_idx in sim.active_paths)
_evolve_branched_ast_operators(sim::BranchedSimulatorOperators, i::Vector{<:Number}) = Dict(path_idx => i for path_idx in sim.active_paths)

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
				if child_head == :break
					# For break statements, clear all active paths to stop the simulation
					# from continuing through the loop
					empty!(sim.active_paths)
				end
				return
			end
			return_value = _evolve_branched_ast_operators(sim, child_expr)
			if child_head == :return
				return return_value
			end
			isempty(sim.active_paths) && return  # No active paths left
		end

		# Subroutines return up to one value of classical type, signified by the return keyword. 
		# If there is no return value, the empty return keyword may be used to immediately exit from the subroutine, which implicitly returns the void type.
	elseif expr_type == :return
		if length(expr.args) == 0 # No return value
			for path_idx in sim.active_paths
				sim.return_values[path_idx] = nothing
			end
		else
			sim.return_values = merge(sim.return_values, _evolve_branched_ast_operators(sim, expr.args[1]))
		end
		sim.active_paths = []

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
			set_variable!(sim, path_idx, var_name, var)
		end

	elseif expr_type ∈ (:continue, :break)
		# Loop control statements
		if expr_type == :break
			# For break statements, clear all active paths to stop the simulation
			# from continuing through the loop
			empty!(sim.active_paths)
		end
		return

	elseif expr_type == :for
		_handle_for_loop(sim, expr)

	elseif expr_type == :switch
		_handle_switch_statement(sim, expr)

	elseif expr_type == :alias
		alias_name = Quasar.name(expr)
		right_hand_side = expr.args[1]

		# Check if the right-hand side is a direct expression or wrapped in args
		if length(right_hand_side.args) > 0
			right_hand_side = right_hand_side.args[1]
			if length(right_hand_side.args) > 0
				right_hand_side = right_hand_side.args[end]
			end
		end

		# Process the alias based on the right-hand side type
		if head(right_hand_side) == :binary_op
			# Handle concatenation
			right_hand_side.args[1] == Symbol("++") || error("Right hand side of alias must be either an identifier, indexed identifier, or concatenation")
			concat_left = right_hand_side.args[2]
			concat_right = right_hand_side.args[3]

			# Check if either side is a physical qubit
			if head(concat_left) == :hw_qubit || head(concat_right) == :hw_qubit
				error("Cannot alias physical qubits. The let keyword allows declared quantum bits and registers to be referred to by another name.")
			end

			# Get qubits from both sides
			left_qs_by_path = evaluate_qubits(sim, concat_left)
			right_qs_by_path = evaluate_qubits(sim, concat_right)

			# Process for each active path
			for path_idx in sim.active_paths
				# Skip if this path doesn't have qubit indices for either side
				if !haskey(left_qs_by_path, path_idx) || !haskey(right_qs_by_path, path_idx)
					continue
				end

				# Concatenate the qubit indices
				left_qs = sort(left_qs_by_path[path_idx])
				right_qs = sort(right_qs_by_path[path_idx])
				alias_qubits = vcat(left_qs, right_qs)

				# Store the alias in variables
				set_variable!(sim, path_idx, alias_name, ClassicalVariable(alias_name, :qubit_alias, alias_qubits, false))
			end

		elseif head(right_hand_side) == :identifier
			referent_name = Quasar.name(right_hand_side)

			# Check if it's a physical qubit
			if head(right_hand_side) == :hw_qubit
				error("Cannot alias physical qubits. The let keyword allows declared quantum bits and registers to be referred to by another name.")
			end

			# Get the qubits from the identifier
			qubits_by_path = evaluate_qubits(sim, right_hand_side)

			# Process for each active path
			for path_idx in sim.active_paths
				# Skip if this path doesn't have qubit indices
				if !haskey(qubits_by_path, path_idx)
					continue
				end

				alias_qubits = qubits_by_path[path_idx]

				# Store the alias in variables
				set_variable!(sim, path_idx, alias_name, ClassicalVariable(alias_name, :qubit_alias, alias_qubits, false))
			end
		elseif head(right_hand_side) == :indexed_identifier
			referent_name = Quasar.name(right_hand_side)

			# Get the qubits from the indexed identifier
			qubits_by_path = evaluate_qubits(sim, right_hand_side)

			# Process for each active path
			for path_idx in sim.active_paths
				# Skip if this path doesn't have qubit indices
				if !haskey(qubits_by_path, path_idx)
					continue
				end

				alias_qubits = qubits_by_path[path_idx]

				# Store the alias in variables
				set_variable!(sim, path_idx, alias_name, ClassicalVariable(alias_name, :qubit_alias, alias_qubits, false))
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

		# Create a dictionary to store results for each path
		results = Dict{Int, Any}()

		# Process each active path
		for path_idx in sim.active_paths
			# First check path-specific variables
			var = get_variable(sim, path_idx, id_name)

			if !isnothing(var)
				# This is because we store the value of singular bits as a bit register with one bit
				if var.type == Quasar.SizedBitVector && length(var.val) == 1
					results[path_idx] = var.val[1]
				elseif !isnothing(var.val)
					results[path_idx] = var.val
				end
				# Then check qubit_mapping
			elseif haskey(sim.qubit_mapping, id_name)
				results[path_idx] = sim.qubit_mapping[id_name]
			else
				# If we get here, the identifier is not defined
				error("No identifier $id_name defined for path $path_idx")
			end
		end

		return results

		# Indexed identifier reference
	elseif expr_type == :indexed_identifier
		identifier_name = Quasar.name(expr)

		# Only process for active paths
		if isempty(sim.active_paths)
			error("No active paths when evaluating indexed identifier $identifier_name")
		end

		indices = _evolve_branched_ast_operators(sim, expr.args[2])


		# Create a dictionary to store results for each path
		results = Dict{Int, Any}()

		# Process each active path
		for (path_idx, index) in indices
			# Get the variable
			var = get_variable(sim, path_idx, identifier_name)

			if !isnothing(var)
				# Check if it's a bit array (SizedBitVector)
				if var.type == Quasar.SizedBitVector
					# Convert to 1-based indexing for Julia arrays
					julia_index = index + 1

					if julia_index <= length(var.val)
						# Access the bit directly from the boolean array
						results[path_idx] = var.val[julia_index]
					end
				elseif !isnothing(var.type) && var.type <: Quasar.SizedArray
					# For other array types
					# Convert to 1-based indexing for Julia arrays
					flat_ix = (index isa Vector) ? index : [index]
					flat_ix = flat_ix .+ 1

					results[path_idx] = var.val[only(flat_ix)]
				end
			else
				# Then check qubit_mapping
				# Temporarily set this path as the only active one for qubit evaluation
				original_active_paths = copy(sim.active_paths)
				sim.active_paths = [path_idx]
				path_results = evaluate_qubits(sim, expr)
				sim.active_paths = original_active_paths

				if haskey(path_results, path_idx)
					results[path_idx] = path_results[path_idx]
				end
			end
		end

		return results

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
		# Return the literal value as a dictionary mapping path indices to the value
		return Dict(path_idx => expr.args[1] for path_idx in sim.active_paths)

	elseif expr_type == :array_literal
		# Evaluate array elements
		results = Dict{Int, Vector{Any}}()

		for element in expr.args
			result_elem = _evolve_branched_ast_operators(sim, element)

			for (path_idx, result) in result_elem
				if path_idx in keys(results)
					push!(results[path_idx], result)
				else
					results[path_idx] = [result]
				end
			end
		end

		return results

	elseif expr_type == :range
		# Evaluate range parameters
		start = _evolve_branched_ast_operators(sim, expr.args[1])
		step = length(expr.args) > 2 ? _evolve_branched_ast_operators(sim, expr.args[2]) : 1
		stop = length(expr.args) > 2 ? _evolve_branched_ast_operators(sim, expr.args[3]) : _evolve_branched_ast_operators(sim, expr.args[2])
		results = Dict{Int, Any}()
		for path_idx in sim.active_paths
			results[path_idx] = collect(start[path_idx]:step[path_idx]:stop[path_idx])
		end
		return results

	elseif expr_type == :binary_op
		# Binary operation
		op = expr.args[1]
		lhs = _evolve_branched_ast_operators(sim, expr.args[2])
		rhs = _evolve_branched_ast_operators(sim, expr.args[3])

		results = Dict{Int, Any}()

		# Get all unique path indices from both operands
		path_indices = keys(lhs)

		# For each path, evaluate the binary operation with path-specific values
		for path_idx in path_indices
			# Get the values for this path
			lhs_val = get(lhs, path_idx, lhs)
			rhs_val = get(rhs, path_idx, rhs)

			# Evaluate the binary operation for this path
			results[path_idx] = evaluate_binary_op(op, lhs_val, rhs_val)
		end

		return results

	elseif expr_type == :unary_op
		# Unary operation
		op = expr.args[1]
		arg = _evolve_branched_ast_operators(sim, expr.args[2])

		results = Dict{Int, Any}()

		# For each path, evaluate the unary operation with path-specific value
		for (path_idx, path_val) in arg
			results[path_idx] = evaluate_unary_op(op, path_val)
		end

		return results


	elseif expr_type == :cast
		casting_to = expr.args[1].args[1]
		value = _evolve_branched_ast_operators(sim, expr.args[2])
		results = Dict{Int, Any}()

		# For each path, apply the cast operation to the path-specific value
		for (path_idx, path_val) in value
			# Apply the appropriate cast based on the target type
			if casting_to == Bool
				results[path_idx] = path_val > 0
			elseif casting_to isa Quasar.SizedInt
				if typeof(path_val) <: AbstractFloat
					results[path_idx] = floor(path_val)
				elseif path_val isa Vector
					results[path_idx] = parse(Int, join(path_val), base = 2)
				else
					results[path_idx] = Int(path_val)
				end
			elseif casting_to isa Quasar.SizedFloat
				results[path_idx] = Float64(path_val)
			elseif casting_to isa Quasar.SizedBitVector
				if typeof(path_val) <: AbstractFloat
					error("Unsupported conversion from float to bit vector")
				elseif typeof(path_val) <: Int
					results[path_idx] = digits(path_val, base = 2)
				else
					error("Unsupported conversion to bit vector from $(typeof(path_val))")
				end
			else
				error("Unsupported cast type: $casting_to")
			end
		end

		return results

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
	qubit_indices_dict = evaluate_qubits(sim, expr.args[1])

	new_active_paths = Int[]

	# Structure to store measurement outcomes for each path
	path_outcomes = Dict{Int, Dict{String, Int}}()

	for path_idx in sim.active_paths
		# Skip if this path doesn't have qubit indices
		if !haskey(qubit_indices_dict, path_idx)
			continue
		end

		qubit_indices = qubit_indices_dict[path_idx]

		# Everytime a measurement occurs, the resulting path indices are added to the following list
		paths_to_measure = [path_idx]

		path_outcomes[path_idx] = Dict{String, Int}()

		# For every active path, measure each qubit in the qubit indices
		for qubit_idx in qubit_indices
			# Find the qubit name for this index
			qubit_name = ""
			for (name, idx) in sim.qubit_mapping
				if idx == qubit_idx
					qubit_name = name
					break
				end
			end

			# If we couldn't find the qubit name directly, check if it's a function parameter
			if qubit_name == ""
				# Check all variables in this path to see if any are qubit parameters
				for (var_name, var) in sim.variables[path_idx]
					if var.type == :qubit_declaration && !isnothing(var.val)
						stored_qubit_name = var.val
						if haskey(sim.qubit_mapping, stored_qubit_name) && sim.qubit_mapping[stored_qubit_name] == qubit_idx
							qubit_name = stored_qubit_name
							break
						end
					end
				end
			end

			# If we still couldn't find the qubit name, generate a default name
			if qubit_name == ""
				error("Qubit name can't be found!")
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
	_deep_copy_variables(variables::Dict{String, FramedVariable})

Make a deep copy of a variables dictionary to avoid aliasing issues with mutable objects.
"""
function _deep_copy_variables(variables::Dict{String, FramedVariable})
	result = Dict{String, FramedVariable}()
	for (name, var) in variables
		# Make a deep copy of the value if it's a collection
		val_copy = var.val isa AbstractArray ? copy(var.val) : var.val
		result[name] = FramedVariable(var.name, var.type, val_copy, var.is_const, var.frame_number)
	end
	return result
end

"""
	_create_block_scope(sim::BranchedSimulatorOperators)

Helper function to create a new scope for block statements (for loops, if/else, while loops).
Unlike function and gate scopes, block scopes inherit all variables from the containing scope.
Returns a dictionary mapping path indices to their original variable dictionaries.
Increments the current frame number to indicate entering a new scope.
"""
function _create_block_scope(sim::BranchedSimulatorOperators)
	original_variables = Dict{Int, Dict{String, FramedVariable}}()

	# Increment the current frame as we're entering a new scope
	sim.curr_frame += 1

	# Save current variables state for all active paths (dont deep copy to include aliasing)
	for path_idx in sim.active_paths
		original_variables[path_idx] = copy(sim.variables[path_idx])
	end

	return original_variables
end

"""
	_handle_for_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a for loop in the branched simulation model.
Creates a new scope for the loop body and properly handles variable scoping.
"""
function _handle_for_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)
	for_loop = convert(Vector{QasmExpression}, expr.args)
	loop_variable_type = for_loop[1].args[1]
	loop_variable_name = for_loop[2].args[1]::String

	# Evaluate loop range
	loop_variable_values_dict = _evolve_branched_ast_operators(sim, for_loop[3])
	loop_body = for_loop[4]::QasmExpression

	# Calculate paths not to add
	paths_not_to_add = setdiff(range(1, length(sim.instruction_sequences)), sim.active_paths)

	# Create a new scope for the loop
	original_variables = _create_block_scope(sim)

	# For each value in the range
	for (path_idx, loop_variable_vals) in loop_variable_values_dict
		# Set active paths to only those that have this loop value
		sim.active_paths = [path_idx]

		for loop_value in loop_variable_vals
			# Set the loop variable for each active path
			for path_idx in sim.active_paths
				set_variable!(sim, path_idx, loop_variable_name, ClassicalVariable(loop_variable_name, loop_variable_type, loop_value, false))
			end

			# Process loop body
			_evolve_branched_ast_operators(sim, loop_body)

			# If no active paths left, continue to next value
			isempty(sim.active_paths) && continue
		end
	end

	# Restore active paths
	sim.active_paths = setdiff(range(1, length(sim.instruction_sequences)), paths_not_to_add)

	# Restore original scope
	_restore_original_scope(sim, original_variables)
end

"""
	_handle_conditional(sim::BranchedSimulatorOperators, expr::QasmExpression)
Handle a conditional statement, filtering paths based on the condition.
Creates a new scope for both the if and else branches.
"""
function _handle_conditional(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# Find if there's an else branch
	has_else = findfirst(e->head(e) == :else, convert(Vector{QasmExpression}, expr.args))
	last_expr = !isnothing(has_else) ? length(expr.args) - 1 : length(expr.args)

	new_paths = []

	# Evaluate condition for each active path
	true_paths = Int[]
	false_paths = Int[]

	condition_values = _evaluate_condition(sim, expr.args[1])

	for (path_idx, condition_value) in condition_values
		if condition_value || condition_value > 0
			push!(true_paths, path_idx)
		else
			push!(false_paths, path_idx)
		end
	end

	# Process if branch for true paths
	if !isempty(true_paths)
		sim.active_paths = true_paths

		# Create a new scope for the if branch
		original_variables = _create_block_scope(sim)

		# Process if branch
		for child_expr in expr.args[2:last_expr]
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && break
		end

		# Restore original scope
		_restore_original_scope(sim, original_variables)

		# Add surviving paths to new_paths
		append!(new_paths, sim.active_paths)
	end

	# Process else branch for false paths
	if !isnothing(has_else) && !isempty(false_paths)
		sim.active_paths = false_paths

		# Create a new scope for the else branch
		original_variables = _create_block_scope(sim)

		# Process else branch
		for child_expr in expr.args[has_else].args
			_evolve_branched_ast_operators(sim, child_expr)
			isempty(sim.active_paths) && break
		end

		# Restore original scope
		_restore_original_scope(sim, original_variables)

		# Add surviving paths to new_paths
		append!(new_paths, sim.active_paths)
	else
		append!(new_paths, false_paths) # Add false paths directly if no else branch
	end

	# Update active paths
	sim.active_paths = new_paths
end


"""
	_evaluate_condition(sim::BranchedSimulatorOperators, path_idx::Int, condition_expr::QasmExpression)

Evaluate a condition expression in the context of a specific path.
This handles path-specific variables and measurement results.
"""
function _evaluate_condition(sim::BranchedSimulatorOperators, condition_expr::QasmExpression)
	# Save the current active paths
	original_active_paths = copy(sim.active_paths)
	result = _evolve_branched_ast_operators(sim, condition_expr)

	# Restore active paths
	sim.active_paths = original_active_paths

	return result
end

"""
	_handle_switch_statement(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a switch statement in the branched simulation model. For more information, visit
https://openqasm.com/language/classical.html#the-switch-statement
"""
function _handle_switch_statement(sim::BranchedSimulatorOperators, expr::QasmExpression)
	all_cases = convert(Vector{QasmExpression}, expr.args[2:end])
	default_idx = findfirst(expr->head(expr) == :default, all_cases)

	# Group paths by their case values
	case_paths = Dict{Any, Vector{Int}}()
	default_paths = Int[]

	# For each active path, evaluate the condition and determine which case it should follow
	original_active_paths = copy(sim.active_paths)

	# Evaluate switch condition for all paths at once
	sim.active_paths = original_active_paths
	case_vals = _evolve_branched_ast_operators(sim, expr.args[1])

	all_cases_vals = []
	# Get all of the case values themselves
	for (case_idx, case) in enumerate(all_cases)
		if head(case) == :case
			push!(all_cases_vals, _evolve_branched_ast_operators(sim, case.args[1]))
		end
	end


	for path_idx in original_active_paths
		# Get the case value for this path
		case_val = get(case_vals, path_idx, nothing)

		# Skip paths with no case value
		isnothing(case_val) && continue

		# Find matching case for this path
		path_matched = false

		for (case_idx, case) in enumerate(all_cases_vals)
			# Case value matched
			if case[path_idx] == case_val
				if !haskey(case_paths, case_idx)
					case_paths[case_idx] = Int[]
				end
				push!(case_paths[case_idx], path_idx)
				path_matched = true
				break
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
Creates a new scope for the loop body and properly handles variable scoping.
"""
function _handle_while_loop(sim::BranchedSimulatorOperators, expr::QasmExpression)
	loop_body = expr.args[2]

	# Calculate paths not to add
	paths_not_to_add = setdiff(range(1, length(sim.instruction_sequences)), sim.active_paths)

	# Create a new scope for the entire while loop
	original_variables = _create_block_scope(sim)

	# Keep track of paths that should continue looping
	continue_paths = copy(sim.active_paths)

	while !isempty(continue_paths)
		# Set active paths to only those that should continue looping
		sim.active_paths = continue_paths

		# Evaluate condition for all paths at once
		condition_results = _evolve_branched_ast_operators(sim, expr.args[1])

		# Determine which paths should continue looping
		new_continue_paths = Int[]

		for path_idx in continue_paths
			if haskey(condition_results, path_idx)
				result = condition_results[path_idx]
				condition_value = result isa Bool ? result : result > 0

				if condition_value
					push!(new_continue_paths, path_idx)
				end
			end
		end

		# If no paths should continue, break
		isempty(new_continue_paths) && break

		# Execute the loop body
		sim.active_paths = new_continue_paths
		_evolve_branched_ast_operators(sim, loop_body)

		# Update continue_paths for next iteration
		continue_paths = copy(sim.active_paths)
	end

	# Restore paths that didn't enter the loop
	sim.active_paths = setdiff(range(1, length(sim.instruction_sequences)), paths_not_to_add)

	# Restore original scope
	_restore_original_scope(sim, original_variables)
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

	# Evaluate the LHS and RHS for all active paths (for assignment operation, we ignore the lhs value)
	lhs_value = op == Symbol("=") ? nothing : _evolve_branched_ast_operators(sim, lhs)
	rhs_value = _evolve_branched_ast_operators(sim, rhs)

	# Check if the right-hand side is a measurement result
	is_measurement = rhs_value isa MeasurementValues

	# Get all paths that were created during RHS evaluation
	current_paths = copy(sim.active_paths)

	# This is for assigning to a bit register or any classical variable
	if head(lhs) == :identifier
		for current_path_idx in current_paths
			lhs_val = !isnothing(lhs_value) && haskey(lhs_value, current_path_idx) ? lhs_value[current_path_idx] : nothing # Default to nothing because if no value then we are declaring the variable

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
						# We're just updating the array with the measurement values
						bit_array = var.val
						if bit_array isa Vector && length(bit_array) == length(bit_values)
							# Update the array directly
							var.val = bit_values
						elseif bit_array isa Int && length(bit_values) == 1
							# For when you are just assigning qubit to a bit, not a register
							var.val = bit_values[1]
						else
							error("Assignment must be to a register of the right size")
						end
					end
				end
			else
				# Standard variable assignment
				# Get the existing variable
				var = get(sim.variables[current_path_idx], var_name, nothing)
				new_val = evaluate_binary_op(op, lhs_val, rhs_value[current_path_idx])

				if !isnothing(var)
					var.val = new_val
				else
					error("Variable doesn't exist")
				end
			end
		end
		# This is for assigning to a specific bit in a register
	elseif head(lhs) == :indexed_identifier
		# Handle indexed assignment
		# Get the index for this path
		index_value = _evolve_branched_ast_operators(sim, lhs.args[2])

		for current_path_idx in current_paths
			lhs_val = !isnothing(lhs_value) && haskey(lhs_value, current_path_idx) ? lhs_value[current_path_idx] : nothing # Default to nothing because if no value then we are declaring the variable

			index = get(index_value, current_path_idx, 0)

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
						if bit_array isa Vector && index + 1 <= length(bit_array)
							# Convert to 1-based indexing for Julia arrays
							julia_index = index + 1
							bit_array[julia_index] = evaluate_binary_op(op, lhs_val, outcome)
						end
					end
				end
			else
				# Standard indexed assignment
				# Get the register variable
				register_var = get(sim.variables[current_path_idx], var_name, nothing)
				new_val = evaluate_binary_op(op, lhs_val, rhs_value[current_path_idx])

				if !isnothing(register_var)
					bit_array = register_var.val
					if bit_array isa Vector && index + 1 <= length(bit_array)
						# Convert to 1-based indexing for Julia arrays
						julia_index = index + 1
						bit_array[julia_index] = new_val
					end
				end
			end
		end
	end
	# Add all the active paths generated in the new_active_paths variable
	append!(new_active_paths, sim.active_paths)
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
			register_size_dict = _evolve_branched_ast_operators(sim, var_type.size)

			# Process each path individually
			for (path_idx, register_size) in register_size_dict
				# If the size is -1, then it is a single bit, not a register
				if register_size == -1
					# Initialize with a single boolean value
					set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, typeof(var_type), [1], false))
				else
					# Initialize the register with a boolean array
					set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, typeof(var_type), [0 for _ in 1:register_size], false))
				end
			end
		else
			# Initialize variable for each active path
			for path_idx in sim.active_paths
				set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, typeof(var_type), nothing, false))
			end
		end
		# This for all other classical variables
	elseif head(expr.args[2]) == :classical_assignment
		var_name = Quasar.name(expr.args[2].args[1].args[2])
		var_type = typeof(expr.args[1].args[1])

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
			set_variable!(sim, path_idx, var_name, ClassicalVariable(var_name, var_type, undef, true))
		end

		# Handle the assignment
		_handle_classical_assignment(sim, expr.args[2])
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
	qubit_name = Quasar.name(expr)
	qubit_size_dict = _evolve_branched_ast_operators(sim, expr.args[2])

	# Get a representative qubit size (should be the same for all paths)
	# Handle any dictionary type, not just Dict{Int, Any}
	qubit_size = if isa(qubit_size_dict, Dict)
		first(values(qubit_size_dict))
	else
		qubit_size_dict
	end

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
	# Save original active paths
	original_active_paths = copy(sim.active_paths)

	# Get modifiers for each path
	modifiers_by_path = evaluate_modifiers(sim, expr)

	# Process each path separately
	for (path_idx, (mod_expr, inner_expr)) in modifiers_by_path
		# Set active paths to just this path
		sim.active_paths = [path_idx]

		# Process gate modifiers for this path
		mods = QasmExpression(:modifiers)
		push!(mods, mod_expr)

		# Get the inner expression
		inner = inner_expr

		# Process nested modifiers
		while head(inner) != :gate_call
			# Get modifiers for this inner expression
			inner_modifiers = evaluate_modifiers(sim, inner)

			# Since we're only processing one path, there should be only one entry
			if haskey(inner_modifiers, path_idx)
				mod_expr, inner = inner_modifiers[path_idx]
				push!(mods, mod_expr)
			else
				# If this path doesn't have modifiers for this inner expression, break
				break
			end
		end

		# Add modifiers to the inner gate call
		push!(inner, mods)

		# Process the gate call with modifiers
		_evolve_branched_ast_operators(sim, inner)
	end

	# Restore original active paths
	sim.active_paths = original_active_paths
end

"""
	_handle_gate_call(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a gate call in the branched simulation model.
Instead of applying gates directly, store them in the instruction sequence.
When a measurement is needed, all stored instructions will be applied.

According to the scoping rules:
- Inside the gate, only const variables and previously defined gates and subroutines are visible
- Variables defined in the gate scope are local to the gate body
- Parameters behave as if they were defined in the gate scope
- The gate's local variables end at the end of the gate body
- The qubit identifiers in a gate definition are valid only within the definition of the gate
"""
function _handle_gate_call(sim::BranchedSimulatorOperators, expr::QasmExpression)
	gate_name = Quasar.name(expr)

	# Extract gate parameters and evaluate them in the current scope
	evaluated_params = Dict{Int, Vector{Any}}()
	if length(expr.args[2].args) > 0 && !isnothing(expr.args[2])
		param_exprs = convert(Vector{QasmExpression}, expr.args[2].args)[1]
		# Evaluate parameters in the current scope - these will be dictionaries mapping path indices to values
		param_results = _evolve_branched_ast_operators(sim, param_exprs)

		# Convert the parameter results to a dictionary mapping path indices to parameter vectors
		for path_idx in sim.active_paths
			if param_results isa Dict && haskey(param_results, path_idx)
				if param_results[path_idx] isa Vector
					evaluated_params[path_idx] = param_results[path_idx]
				else
					evaluated_params[path_idx] = [param_results[path_idx]]
				end
			end
		end
	end

	# Extract target qubit indices and evaluate them in the current scope
	target_indices_dict = Dict{Int, Vector{Int}}()
	if length(expr.args) > 2 && !isnothing(expr.args[3])
		if head(expr.args[3]) == :qubit_targets
			for target in expr.args[3].args
				qubit_indices = evaluate_qubits(sim, target)
				for (path_idx, indices) in qubit_indices
					if !haskey(target_indices_dict, path_idx)
						target_indices_dict[path_idx] = Int[]
					end
					append!(target_indices_dict[path_idx], indices)
				end
			end
		else
			# Handle the case where it's not a :qubit_targets expression
			qubit_indices = evaluate_qubits(sim, expr.args[3])
			for (path_idx, indices) in qubit_indices
				target_indices_dict[path_idx] = indices
			end
		end
	end

	# Check if there are any modifiers to apply
	modifiers = Vector{QasmExpression}()
	if length(expr.args) > 3 && !isnothing(expr.args[4]) && head(expr.args[4]) == :modifiers
		modifiers = expr.args[4].args
	end

	# Process each active path
	for path_idx in sim.active_paths
		# Skip if this path doesn't have target indices
		if !haskey(target_indices_dict, path_idx)
			continue
		end

		# Get target indices for this path
		target_indices = target_indices_dict[path_idx]

		# Get gate operator
		if haskey(sim.gate_defs, gate_name)
			gate_def = sim.gate_defs[gate_name]

			# Check if it's a built-in gate (in the mapping) or a user-defined gate
			if haskey(sim.gate_name_mapping, gate_name)
				# Built-in gate
				path_params = haskey(evaluated_params, path_idx) ? evaluated_params[path_idx] : []

				# Special case for GPhase gate which needs to know the number of qubits
				if gate_name == "gphase"
					# Get the number of target qubits
					N = sim.n_qubits-length(target_indices)

					# Create GPhase{N} gate with the correct number of qubits
					gate_op = if !isempty(path_params)
						GPhase{N}(tuple(path_params...))
					else
						GPhase{N}()
					end
				else
					# Built-in gate - use the mapping to get the gate type
					gate_type = getfield(BraketSimulator, sim.gate_name_mapping[gate_name])

					# Create gate operator
					gate_op = if !isempty(path_params)
						gate_type(tuple(path_params...))
					else
						gate_type()
					end
				end

				# Apply modifiers if any
				if !isempty(modifiers)
					gate_op = _apply_modifiers(gate_op, modifiers, target_indices)
				end

				# Store gate instruction for this path
				instruction = Instruction(gate_op, target_indices)
				push!(sim.instruction_sequences[path_idx], instruction)
			else
				# User-defined gate - create a new scope and execute the body
				# Save original active paths
				original_active_paths = copy(sim.active_paths)

				# Set active paths to just this path
				sim.active_paths = [path_idx]

				# Create a new scope with only const variables
				original_variables = _create_const_only_scope(sim)

				# Bind parameters to the gate scope if this is a parametric gate
				if haskey(evaluated_params, path_idx) && !isempty(gate_def.arguments)
					path_params = evaluated_params[path_idx]
					for (i, param_name) in enumerate(gate_def.arguments)
						if i <= length(path_params)
							# Store parameters with their actual type, not just Any
							param_value = path_params[i]
							param_type = typeof(param_value)
							set_variable!(sim, path_idx, param_name, param_type, param_value, false)
						end
					end
				end

				# Bind qubit targets to the gate scope
				if !isempty(target_indices) && !isempty(gate_def.qubit_targets)
					for (i, target_name) in enumerate(gate_def.qubit_targets)
						if i <= length(target_indices)
							set_variable!(sim, path_idx, target_name, :qubit_declaration, target_indices[i], false)
						end
					end
				end

				# Execute the gate body
				instruction = gate_def.body
				if instruction isa Quasar.QasmExpression
					for gate_call in instruction.args
						_evolve_branched_ast_operators(sim, gate_call)
					end
				else
					# For non-QasmExpression instructions (e.g., direct gate operators)
					gate_op = StructTypes.constructfrom(QuantumOperator, instruction)

					# Apply modifiers if any
					if !isempty(modifiers)
						gate_op = _apply_modifiers(gate_op, modifiers, target_indices)
					end

					# Store gate instruction for this path
					gate_instruction = Instruction(gate_op, target_indices)
					push!(sim.instruction_sequences[path_idx], gate_instruction)
				end

				# Restore original variables
				_restore_original_scope(sim, original_variables)

				# Restore original active paths
				sim.active_paths = original_active_paths
			end
		else
			error("Gate $gate_name not defined!")
		end
	end

	return nothing
end

"""
	_handle_gate_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a gate definition in the branched simulation model.

The definitions of gates (gate) introduce a new scope. The gate statements are only valid directly 
within the global scope of the program.

Inside the definition of the gate, symbols that were already defined in the global scope with the const modifier, 
or previously defined gates and subroutines are visible. Globally scoped variables without the const modifier 
are not visible inside the definition. In other words, gates cannot close over variables that may be modified at run-time.

Variables defined in the parameter specifications of gates behave for scoping purposes as if they were defined 
in the scope of the definition. The lifetime of these local variables ends at the end of the gate body, and they 
are not accessible after the gate body. Similarly, the qubit identifiers in a gate definition are valid only 
within the definition of the gate.

The identifier of a gate is available in the scope of its own body, allowing direct recursion. For gates, 
the direct recursion is unlikely to ever be useful, since this would generally be non-terminating.

Local gate variables, including parameters and qubit definitions, may shadow variables defined in the outer scope. 
Inside the body, the identifier will refer to the local variable instead. After the definition of the body has 
completed (and we are back in the global scope), the identifier will refer to the same variable it did before the gate.

Aliases can be declared within gate scopes, and have the same lifetime and visibility as other local variables.
"""
function _handle_gate_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)
	# Process gate definition following the provided pattern
	gate_def         = expr.args
	gate_name        = Quasar.name(expr)
	gate_arguments   = gate_def[2]::QasmExpression
	gate_def_targets = gate_def[3]::QasmExpression
	gate_body        = gate_def[4]::QasmExpression

	# Extract argument expressions and names
	argument_exprs = !isempty(gate_arguments.args) ? convert(Vector{QasmExpression}, gate_arguments.args[1]) : QasmExpression[]
	argument_names = String[]
	for arg in argument_exprs
		if arg.args[1] isa String
			push!(argument_names, arg.args[1])
		else
			push!(argument_names, string(arg.args[1]))
		end
	end

	# Extract qubit targets
	qubit_targets = String[]
	if !isempty(gate_def_targets.args)
		for target in convert(Vector{QasmExpression}, gate_def_targets.args[1])
			push!(qubit_targets, Quasar.name(target))
		end
	end

	# Store the gate definition in the simulator
	sim.gate_defs[gate_name] = GateDefinition(gate_name, argument_names, qubit_targets, gate_body)
end

"""
	_create_const_only_scope(sim::BranchedSimulatorOperators)

Helper function to create a new scope that only includes const variables from the current scope.
Returns a dictionary mapping path indices to their original variable dictionaries.
Increments the current frame number to indicate entering a new scope.
"""
function _create_const_only_scope(sim::BranchedSimulatorOperators)
	original_variables = Dict{Int, Dict{String, FramedVariable}}()

	# Increment the current frame as we're entering a new scope
	sim.curr_frame += 1

	# Current frame that variables in this scope will be assigned to
	current_frame = sim.curr_frame

	# Save current variables state and create new scopes with only const variables
	for path_idx in sim.active_paths
		original_variables[path_idx] = copy(sim.variables[path_idx])

		# Create a new variable scope
		new_scope = Dict{String, FramedVariable}()

		# Copy only const variables to the new scope
		for (var_name, var) in sim.variables[path_idx]
			if var.is_const
				new_scope[var_name] = var
			end
		end

		# Update the path's variables to the new scope
		sim.variables[path_idx] = new_scope
	end

	return original_variables
end

"""
	_restore_original_scope(sim::BranchedSimulatorOperators, original_variables::Dict{Int, Dict{String, FramedVariable}})

Helper function to restore the original scope after executing in a temporary scope.
For paths that existed before the function call, restore the original scope with original values.
For new paths created during the function call, remove all variables that were instantiated in the current frame.

This implementation handles scoping for newly generated paths by tracking the frame where each variable was declared.
It also properly handles variable shadowing by restoring the original variables when exiting a scope.
"""
function _restore_original_scope(sim::BranchedSimulatorOperators, original_variables::Dict{Int, Dict{String, FramedVariable}})
	# Get all paths that existed before the function call
	original_paths = collect(keys(original_variables))

	# Store the current frame that we're exiting from
	exiting_frame = sim.curr_frame

	# Decrement the current frame as we're exiting a scope
	sim.curr_frame -= 1

	# For paths that existed before, restore the original scope
	for path_idx in sim.active_paths
		if haskey(original_variables, path_idx)
			# Create a new scope that combines original variables with updated values
			new_scope = Dict{String, FramedVariable}()

			# First, copy all original variables to ensure we don't lose any
			for (var_name, orig_var) in original_variables[path_idx]
				new_scope[var_name] = orig_var
			end

			# Then update any variables that were modified in outer scopes
			for (var_name, current_var) in sim.variables[path_idx]
				if current_var.frame_number < exiting_frame && haskey(new_scope, var_name)
					# This is a variable from an outer scope that was modified
					# Keep the original variable's frame number but use the updated value
					orig_var = new_scope[var_name]
					new_scope[var_name] = FramedVariable(
						orig_var.name,
						orig_var.type,
						copy(current_var.val),  # Use the updated value
						orig_var.is_const,
						orig_var.frame_number,  # Keep the original frame number
					)
				elseif current_var.frame_number < exiting_frame && !haskey(new_scope, var_name)
					# This is a new variable declared in an outer scope
					new_scope[var_name] = current_var
				end
				# Variables declared in the current frame (frame_number == exiting_frame) are discarded
			end

			# Update the path's variables to the new scope
			sim.variables[path_idx] = new_scope
		else
			# This is a new path created during function execution or measurement
			# We need to keep variables from outer scopes but remove variables from the current frame

			# Create a new scope for this path
			new_scope = Dict{String, FramedVariable}()

			# Find a reference path to copy variables from
			reference_path = original_paths[1]

			# Copy all variables from the current path that were declared in outer frames
			for (var_name, var) in sim.variables[path_idx]
				if var.frame_number < exiting_frame
					# This variable was declared in an outer scope, keep it
					new_scope[var_name] = var
				end
			end

			# Also copy variables from the reference path that might not be in this path
			# This ensures that all paths have the same variable names after exiting a scope
			for (var_name, var) in original_variables[reference_path]
				if !haskey(new_scope, var_name)
					# Create a copy of the variable with the same frame number
					new_scope[var_name] = FramedVariable(
						var.name,
						var.type,
						copy(var.val),
						var.is_const,
						var.frame_number,
					)
				end
			end

			# Update the path's variables to the new scope
			sim.variables[path_idx] = new_scope
		end
	end
end

"""
	_evaluate_arguments(sim::BranchedSimulatorOperators, arg_exprs::Vector{QasmExpression})

Helper function to evaluate arguments in the current scope before creating a new scope.
Returns a vector of tuples containing the original expression and its evaluated value.
"""
function _evaluate_arguments(sim::BranchedSimulatorOperators, arg_exprs::Vector{QasmExpression})
	evaluated_args = []

	for arg_expr in arg_exprs
		if head(arg_expr) == :indexed_identifier
			# For indexed identifiers like q[1], evaluate the index
			qubit_name = Quasar.name(arg_expr)
			index = _evolve_branched_ast_operators(sim, arg_expr.args[2])
			push!(evaluated_args, (arg_expr, index))
		else
			# For other arguments, evaluate them directly
			arg_value = _evolve_branched_ast_operators(sim, arg_expr)
			push!(evaluated_args, (arg_expr, arg_value))
		end
	end

	return evaluated_args
end

"""
	_bind_qubit_parameter(sim::BranchedSimulatorOperators, path_idx::Int, param_name::String, arg_expr::QasmExpression, arg_value::Any)

Helper function to bind a qubit parameter to a variable in the current scope.
"""
function _bind_qubit_parameter(sim::BranchedSimulatorOperators, path_idx::Int, param_name::String, arg_expr::QasmExpression, arg_value::Any)
	if head(arg_expr) == :indexed_identifier
		# For indexed identifiers like q[1], use the evaluated index
		qubit_name = Quasar.name(arg_expr)
		index = arg_value
		indexed_name = "$qubit_name[$index]"

		if haskey(sim.qubit_mapping, indexed_name)
			qubit_idx = sim.qubit_mapping[indexed_name]
			# Store the qubit index directly
			sim.variables[path_idx][param_name] = ClassicalVariable(param_name, :qubit_declaration, qubit_idx, false)
		else
			error("Qubit $indexed_name not found in qubit mapping")
		end
	else
		# For direct qubit references
		qubit_name = Quasar.name(arg_expr)
		if haskey(sim.qubit_mapping, qubit_name)
			qubit_idx = sim.qubit_mapping[qubit_name]
			# Store the qubit index directly
			sim.variables[path_idx][param_name] = ClassicalVariable(param_name, :qubit_declaration, qubit_idx, false)
		else
			error("Qubit $qubit_name not found in qubit mapping")
		end
	end
end

"""
	_handle_function_call(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a function call in the branched simulation model.

According to the scoping rules:
- Inside the function, only const variables and previously defined gates and subroutines are visible
- Variables defined in the function scope are local to the function body
- Parameters behave as if they were defined in the function scope
- The function's local variables end at the end of the function body
- Subroutines cannot contain qubit declarations in their bodies
"""
function _handle_function_call(sim::BranchedSimulatorOperators, expr::QasmExpression)
	function_name = Quasar.name(expr)

	# Check if it's a user-defined function
	if haskey(sim.function_defs, function_name)
		# Get the function definition
		func_def = sim.function_defs[function_name]

		# Extract function arguments from the call and evaluate them in the current scope
		concrete_arguments = []
		if length(expr.args) > 1 && !isnothing(expr.args[2])
			concrete_arguments = expr.args[2].args[1]
		end

		# Evaluate all arguments in the current scope - these will be dictionaries mapping path indices to values
		evaluated_arguments = _evolve_branched_ast_operators(sim, concrete_arguments)

		# Create a new scope with only const variables
		original_variables = _create_const_only_scope(sim)

		# Extract parameter definitions from the function definition
		param_defs = func_def.arguments.args[1]
		if head(param_defs) == :array_literal # Deals with only one parameter input
			param_defs = param_defs.args
		else
			param_defs = [param_defs]
		end

		# Bind arguments to parameters if we have parameter definitions
		for path_idx in sim.active_paths
			arguments = evaluated_arguments[path_idx]

			for (i, param_def) in enumerate(param_defs)
				param_name = Quasar.name(param_def)
				param_type = head(param_def)
				arg_expr = expr.args[2].args[1]
				if head(arg_expr) == :array_literal
					arg_expr = arg_expr.args[i]
				end
				path_arg_value = arguments[i]

				if param_type == :qubit_declaration
					# For qubit parameters, we need to use the qubit index
					if head(arg_expr) == :indexed_identifier
						# For indexed identifiers like q[1], use the evaluated index
						qubit_name = Quasar.name(arg_expr)
						index = path_arg_value[1]
						indexed_name = "$qubit_name[$index]"

						if haskey(sim.qubit_mapping, indexed_name)
							qubit_idx = sim.qubit_mapping[indexed_name]
							# Store the qubit index directly
							sim.variables[path_idx][param_name] = FramedVariable(param_name, :qubit_declaration, qubit_idx, false, sim.curr_frame+1)
						else
							error("Qubit $indexed_name not found in qubit mapping")
						end
					else
						# For direct qubit references
						qubit_name = Quasar.name(arg_expr)
						if haskey(sim.qubit_mapping, qubit_name)
							qubit_idx = sim.qubit_mapping[qubit_name]
							# Store the qubit index directly
							sim.variables[path_idx][param_name] = FramedVariable(param_name, :qubit_declaration, qubit_idx, false, sim.curr_frame+1)
						else
							error("Qubit $qubit_name not found in qubit mapping")
						end
					end
				elseif param_type == :classical_declaration
					var_type = param_def.args[1].args[1]

					# For bit parameters, update the value to a vector if not already
					if !(path_arg_value isa Vector) && var_type == Quasar.SizedBitVector
						path_arg_value = [path_arg_value]
					end
					sim.variables[path_idx][param_name] = FramedVariable(param_name, param_type, path_arg_value, false, sim.curr_frame+1)
				else
					# For other parameter types, use the evaluated value
					sim.variables[path_idx][param_name] = FramedVariable(param_name, param_type, path_arg_value, false, sim.curr_frame+1)
				end
			end
		end

		# Make sure we're handling the function body correctly
		if func_def.body isa Vector{QasmExpression}
			# If it's a vector of expressions, process each one
			for expr in func_def.body
				_evolve_branched_ast_operators(sim, expr)
			end
		else
			# Otherwise, process it directly
			_evolve_branched_ast_operators(sim, func_def.body)
		end

		active_paths = keys(sim.return_values)
		append!(sim.active_paths, active_paths)

		# Restore original variables
		_restore_original_scope(sim, original_variables)

		# Return the function's return value
		return sim.return_values
	else
		error("Function $function_name was never defined")
	end
end

"""
	_handle_function_definition(sim::BranchedSimulatorOperators, expr::QasmExpression)

Handle a function definition in the branched simulation model.

The definitions of subroutines (def) introduce a new scope. The def statements are only valid directly 
within the global scope of the program.

Inside the definition of the subroutine, symbols that were already defined in the global scope with the const modifier, 
or previously defined gates and subroutines are visible. Globally scoped variables without the const modifier 
are not visible inside the definition. In other words, subroutines cannot close over variables that may be modified at run-time.

Variables defined in subroutine scopes are local to the subroutine body. Variables defined in the parameter specifications 
of subroutines behave for scoping purposes as if they were defined in the scope of the definition. The lifetime of these 
local variables ends at the end of the function body, and they are not accessible after the subroutine body.

The identifier of a subroutine is available in the scope of its own body, allowing direct recursion.

Local subroutine variables, including parameters, may shadow variables defined in the outer scope. Inside the body, 
the identifier will refer to the local variable instead. After the definition of the body has completed (and we are 
back in the global scope), the identifier will refer to the same variable it did before the subroutine.

Subroutines cannot contain qubit declarations in their bodies, but can accept variables of type qubit in their parameter lists.
Aliases can be declared within subroutine scopes, and have the same lifetime and visibility as other local variables.
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
