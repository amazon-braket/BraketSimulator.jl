using Quasar: ClassicalVariable, FunctionDefinition, AbstractGateDefinition

# Function to mark a qubit as measured
function mark_measured(qubit::Qubit)
    qubit.measured = true
    return qubit
end

# Function to mark a qubit as unmeasured
function mark_unmeasured(qubit::Qubit)
    qubit.measured = false
    return qubit
end

"""
	BranchedSimulator

A simulator that supports multiple execution paths resulting from measurements.
This implementation stores only instruction sequences, not statevectors.
When a measurement is needed, the state is calculated from scratch by applying
all instructions from the beginning.
"""
mutable struct BranchedSimulator <: AbstractSimulator
	# Base simulator type
	base_simulator::AbstractSimulator

	# Collection of instruction sequences for different paths
	# Each path stores a list of instructions to apply
	states::Vector{Any}
	
	# Active paths (indices into instruction_sequences that are still evolving)
	active_paths::Vector{Int}

	# Tracks and centralizes the measurement outcomes of qubits in MCM
	# Maps path_idx => Dict(qubit_name => measurement_outcome(s))
	measurements::Vector{Dict{String, Union{Int, Vector{Int}}}}

	# Classical variable values for each path
	variables::Vector{Dict{String, ClassicalVariable}}

	# Total number of defined qubits present in each simulation.
	n_qubits::Vector{Int}

	# Number of shots for each simulation
	shots::Vector{Int}

	# Maps qubit names to Qubit objects for each path
	# For registers, stores a vector of Qubit objects
	qubit_mapping::Vector{Dict{String, Union{Qubit, Vector{Qubit}}}}

	# Maps function names to FunctionDefinition objects
	function_defs::Dict{String, FunctionDefinition}

	# Maps gate names to GateDefinition objects
	gate_defs::Dict{String, AbstractGateDefinition}

	# Input parameters for the circuit
	inputs::Dict{String, Any}

	# Constructor
	function BranchedSimulator(simulator::AbstractSimulator; inputs::Dict{String, Any} = Dict{String, Any}())
		# Initialize qubit mapping for each path
		initial_mapping = Dict{String, Union{Qubit, Vector{Qubit}}}()

		# Initialize function and gate definitions
		function_defs = Dict{String, FunctionDefinition}()
		gate_defs = Dict{String, AbstractGateDefinition}()

		# Initialize with built-in gates
		if @isdefined(builtin_gates)
			gate_defs = builtin_gates()
		end

		# Create an empty state of the same type as the simulator's state
		# This ensures type consistency when operations are applied later
		state_type = typeof(get_state(simulator))
		empty_state = similar(get_state(simulator), 0)
		
		# Initialize with empty instruction sequence
		instructions = Vector{Instruction}()

		new(
			simulator,
			[instructions],  # Initialize with empty instruction sequence
			[empty_state],   # Initialize with an empty state of the same type
			[1],
			[Dict{String, Int}()],  # Initialize measurements with qubit name keys
			[Dict{String, ClassicalVariable}()],
			[0],
			[simulator.shots],
			[initial_mapping],       # Initial mapping is identity
			function_defs,
			gate_defs,
			inputs,
		)
	end
end

# Required interface methods
qubit_count(sim::BranchedSimulator) = maximum(sim.n_qubits)
device_id(sim::BranchedSimulator) = device_id(sim.base_simulator)


"""
	update_mapping!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit, outcome::Int)

Updates the qubit mapping structure when a qubit is measured.
Marks the qubit as measured and stores the measurement outcome.
Also updates the indices of the remaining qubits to account for state reduction.
"""
function update_mapping!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit, outcome::Int)
	state = sim.states[path_idx]
	
	# Handle empty state vector case
	if isempty(state) || length(state) == 0
		# Mark the qubit as measured
		mark_measured(qubit)
		
		# Store the measurement outcome in measurements
		sim.measurements[path_idx][qubit.name] = outcome
		
		# Update the qubit count for this path (already 0)
		return
	end
	
	n_qubits = Int(log2(length(state)))

	# If there no qubits present in the state representation, the operation is not allowed,
	# so we will throw an error
	if n_qubits == 0
		throw("There are no qubits to be measured in the state vector")
	end

	# Get the current position of the qubit in the state
	current_qubit_idx = qubit.index
	
	# Mark the qubit as measured
	mark_measured(qubit)

	# Update the qubit mapping
	mapping = sim.qubit_mapping[path_idx]
	mapping[qubit.name] = qubit
	
	# Then, update the indices of the remaining qubits
	# After measuring a qubit, the state is reduced and all higher qubit positions shift down
	for (orig_q, curr_q) in copy(mapping)
		if curr_q isa Qubit && !curr_q.measured && curr_q.index > current_qubit_idx
			curr_q.index -= 1
		elseif curr_q isa Vector{Qubit}
			# Handle qubit registers
			for (i, q) in enumerate(curr_q)
				if !q.measured && q.index > current_qubit_idx
					# This qubit's index needs to be decremented
					curr_q[i].index -= 1
				end
			end
		end
	end

	# Store the measurement outcome in measurements
	sim.measurements[path_idx][qubit.name] = outcome

	# Update the qubit count for this path
	sim.n_qubits[path_idx] = n_qubits - 1
end

# Backward compatibility method
function update_mapping!(sim::BranchedSimulator, path_idx::Int, qubit::Int, outcome::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	update_mapping!(sim, path_idx, qubit_obj, outcome)
end

"""
	get_qubit_index(sim::BranchedSimulator, path_idx::Int, qubit::Qubit) -> Int

Get the current index of a qubit in the state vector.
Handles both individual qubits and qubits that are part of registers.
"""
function get_qubit_index(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	mapping = sim.qubit_mapping[path_idx]

	# Check if this qubit is directly mapped by string key
	qubit_str = qubit.name
	if haskey(mapping, qubit_str)
		value = mapping[qubit_str]
		if value isa Qubit
			return value.index
		elseif value isa Vector{Qubit} && length(value) == 1
			return value[1].index
		end
	end

	# Check if this qubit is part of a register
	for (name, qubits) in mapping
		if qubits isa Vector{Qubit}
			# Try to extract the register base name and check if it's a number
			# This handles both string-based register names and numeric indices
			try
				if occursin("[", name)
					register_base = parse(Int, split(name, "[")[1])
					if register_base <= qubit < register_base + length(qubits)
						return qubits[qubit-register_base+1].index
					end
				end
			catch
				# If parsing fails, it's not a numeric register base, so continue
				continue
			end
		end
	end

	# Default to the qubit index itself if not found in mapping
	return qubit
end

"""
	get_qubit_object(sim::BranchedSimulator, path_idx::Int, qubit::Int) -> Qubit

Get the Qubit object for a given qubit index.
Returns a new Qubit object if the qubit is not found in the mapping.
"""
function get_qubit_object(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	mapping = sim.qubit_mapping[path_idx]

	# Check if this qubit is directly mapped by string key
	qubit_str = string(qubit)
	if haskey(mapping, qubit_str)
		value = mapping[qubit_str]
		if value isa Qubit
			return value
		elseif value isa Vector{Qubit} && length(value) == 1
			return value[1]
		end
	end

	# Check if this qubit is part of a register
	for (name, qubits) in mapping
		if qubits isa Vector{Qubit}
			# Try to extract the register base name and check if it's a number
			# This handles both string-based register names and numeric indices
			try
				if occursin("[", name)
					register_base = parse(Int, split(name, "[")[1])
					if register_base <= qubit < register_base + length(qubits)
						return qubits[qubit-register_base+1]
					end
				end
			catch
				# If parsing fails, it's not a numeric register base, so continue
				continue
			end
		end
	end

	# Create a new Qubit object if not found in mapping
	return Qubit(qubit, string(qubit), false)
end

"""
	expand_state!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)

Expand the state to reincorporate a previously measured qubit.
This is needed when an operation is applied to a qubit that was previously measured.
Uses multiple dispatch to call the appropriate simulator-specific implementation.
"""
function expand_state!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	# Check if the qubit was previously measured
	if !qubit.measured
		return nothing  # Qubit is not measured, nothing to do
	end

	# Get the measurement outcome directly using the qubit name
	outcome = get(sim.measurements[path_idx], qubit.name, -1)
	
	if outcome == -1
		error("Could not find measurement outcome for qubit $(qubit.name)")
	end

	state = sim.states[path_idx]

	# Use multiple dispatch to call the appropriate implementation
	state = expand_state(state, qubit.index, outcome)

	# Mark the qubit as unmeasured
	mark_unmeasured(qubit)

	# Remove the measurement outcome
	delete!(sim.measurements[path_idx], qubit.name)

	# Update the qubit count for this path
	sim.n_qubits[path_idx] = Int(log2(length(state)))
end

# Backward compatibility method
function expand_state!(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	expand_state!(sim, path_idx, qubit_obj)
end

"""
	resize_simulator!(sim::BranchedSimulator, new_qubit_count::Int)

Resize the simulator to accommodate a new total number of qubits.
This is used when dynamically adding qubits during simulation.
"""
function resize_simulator!(sim::BranchedSimulator, path_idx::Int, new_qubit_count::Int)
	# Only resize if we need more qubits
	if all(n_qubits -> new_qubit_count <= n_qubits, sim.n_qubits)
		return nothing
	end

	path_n_qubits = sim.n_qubits[path_idx]

	if new_qubit_count <= path_n_qubits
		return
	end

	# Calculate how many new qubits we need to add for this path
	qubits_to_add = new_qubit_count - path_n_qubits

	state = sim.states[path_idx]

	# Update the state using multiple dispatch
	sim.states[path_idx] = add_qubits(state, qubits_to_add)

	# No need to update the qubit mapping for the new qubits
	# The qubit mappings will be updated when the qubits are actually used

	# Update the qubit count for this path
	sim.n_qubits[path_idx] = new_qubit_count

	# Also update the base simulator's qubit count if possible
	if hasproperty(sim.base_simulator, :qubit_count)
		sim.base_simulator.qubit_count = new_qubit_count
	end
end

"""
	handle_previously_measured_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)

Handle the case where a qubit was previously measured and removed.
Returns true if the qubit was previously measured, false otherwise.
"""
function handle_previously_measured_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	# Debug output
	println("DEBUG: handle_previously_measured_qubit! for qubit $(qubit.name), measured: $(qubit.measured)")
	
	# Check if the qubit is measured
	return qubit.measured
end

# Backward compatibility method
function handle_previously_measured_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	return handle_previously_measured_qubit!(sim, path_idx, qubit_obj)
end

"""
	handle_branching_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit, 
								 probs::Vector{Float64})

Handle the case where measurement can result in multiple outcomes, creating a new branch.
"""
function handle_branching_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit,
	probs::Vector{Float64})
	# Create new path for outcome 1
	original_state = sim.states[path_idx]

	# Calculate shots for each path based on probabilities
	path_shots = sim.shots[path_idx]
	shots_for_outcome_1 = round(Int, path_shots * probs[2])
	shots_for_outcome_0 = path_shots - shots_for_outcome_1

	# Create new state for outcome 1
	new_state = copy(original_state)
	new_state = apply_projection(new_state, qubit.index, 1)

	# Create new path for outcome 1
	new_path_idx = length(sim.states) + 1
	push!(sim.states, new_state)
	push!(sim.shots, shots_for_outcome_1)

	# Copy measurements and variables for the new path
	new_measurements = copy(sim.measurements[path_idx])
	push!(sim.measurements, new_measurements)

	new_variables = copy(sim.variables[path_idx])
	push!(sim.variables, new_variables)

	# Copy qubit mapping and n_qubits for the new path
	push!(sim.qubit_mapping, copy(sim.qubit_mapping[path_idx]))
	push!(sim.n_qubits, sim.n_qubits[path_idx])

	# Create a copy of the qubit for the new path
	new_qubit = Qubit(qubit.index, qubit.name, qubit.measured)

	# Reduce the state for the new path
	update_mapping!(sim, new_path_idx, new_qubit, 1)

	# Update original path for outcome 0
	sim.shots[path_idx] = shots_for_outcome_0
	sim.states[path_idx] = apply_projection(copy(original_state), qubit.index, 0)

	# Reduce the state for the original path
	update_mapping!(sim, path_idx, qubit, 0)
end

# Backward compatibility method
function handle_branching_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Int,
	current_qubit_idx::Int, probs::Vector{Float64})
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	handle_branching_measurement!(sim, path_idx, qubit_obj, probs)
end

"""
	handle_deterministic_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit, 
									 probs::Vector{Float64})

Handle the case where measurement has only one possible outcome.
"""
function handle_deterministic_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit,
	probs::Vector{Float64})
	# Determine the only possible outcome
	outcome = probs[2] > 0.5 ? 1 : 0

	# Apply projection to the state
	sim.states[path_idx] = apply_projection(copy(sim.states[path_idx]), qubit.index, outcome)

	# Update mapping to reflect the measured qubit
	update_mapping!(sim, path_idx, qubit, outcome)
end

# Backward compatibility method
function handle_deterministic_measurement!(sim::BranchedSimulator, path_idx::Int, qubit::Int,
	current_qubit_idx::Int, probs::Vector{Float64})
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	handle_deterministic_measurement!(sim, path_idx, qubit_obj, probs)
end

"""
	measure_qubit(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)

Measure a qubit on a specific path, potentially creating new simulation paths.
Calculates the current state from scratch by applying all instructions.
Returns true if the measurement had two possible outcomes and false otherwise
"""
function measure_qubit(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	# Check if this qubit was previously measured and removed
	if handle_previously_measured_qubit!(sim, path_idx, qubit)
		return false
	end

	# Calculate current state by applying all instructions from scratch
	current_state = calculate_current_state(sim, path_idx)
	
	# Get measurement probabilities for this path
	probs = get_measurement_probabilities(current_state, qubit.index)

	# Determine if we need to branch (both outcomes possible) or not
	if probs[1] > 0 && probs[2] > 0
		# Both outcomes are possible - create a branch
		handle_branching_measurement!(sim, path_idx, qubit, probs)
		return true
	else
		# Only one outcome is possible - no branching needed
		handle_deterministic_measurement!(sim, path_idx, qubit, probs)
		return false
	end
end

# Backward compatibility method
function measure_qubit(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	return measure_qubit(sim, path_idx, qubit_obj)
end

"""
	measure_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)

Measure a qubit on a specific path, creating a temporary new_active_paths vector.
This is a convenience method that creates a new vector for new_active_paths.
"""
function measure_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	# Check if the qubit was previously measured
	if handle_previously_measured_qubit!(sim, path_idx, qubit)
		# If the qubit was already measured, we just need to return the measurement outcome
		return sim.measurements[path_idx][qubit.name]
	end
	
	# If the qubit was not previously measured, measure it now
	path_split = measure_qubit(sim, path_idx, qubit)
	
	# If this is called directly, we need to update the active paths
	if path_split && path_idx in sim.active_paths
		# Remove the original path
		filter!(p -> p != path_idx, sim.active_paths)
		# Add the new paths
		push!(sim.active_paths, path_idx)
		push!(sim.active_paths, length(sim.states))
	end
	
	return nothing
end

# Backward compatibility method
function measure_qubit!(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	return measure_qubit!(sim, path_idx, qubit_obj)
end

"""
	filter_paths!(sim::BranchedSimulator, condition_func)

Filter active paths based on a condition function.
"""
function filter_paths!(sim::BranchedSimulator, condition_func)
	sim.active_paths = filter(path_idx -> condition_func(path_idx, sim), sim.active_paths)
	return nothing
end

"""
	get_measurement_result(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)

Get the measurement result for a specific qubit in a specific path.
"""
function get_measurement_result(sim::BranchedSimulator, path_idx::Int, qubit::Qubit)
	# Return the measurement outcome
	return get(sim.measurements[path_idx], qubit.name, -1)
end

# Backward compatibility method
function get_measurement_result(sim::BranchedSimulator, path_idx::Int, qubit::Int)
	qubit_obj = get_qubit_object(sim, path_idx, qubit)
	return get_measurement_result(sim, path_idx, qubit_obj)
end

"""
	set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, value)

Set a classical variable value for a specific path.
"""
function set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, value::ClassicalVariable)
	sim.variables[path_idx][var_name] = value
end

"""
	get_variable(sim::BranchedSimulator, path_idx::Int, var_name::String, default=nothing)

Get a classical variable value for a specific path.
"""
function get_variable(sim::BranchedSimulator, path_idx::Int, var_name::String, default = nothing)
	return haskey(sim.variables[path_idx], var_name) ? sim.variables[path_idx][var_name].val : default
end

# Helper functions for simulator state manipulation
"""
	get_state(simulator::AbstractSimulator)

Extract the state vector or density matrix from the simulator.
This is a placeholder - actual implementation would depend on simulator type.
"""
function get_state(simulator::AbstractSimulator)
	if hasproperty(simulator, :state_vector)
		return simulator.state_vector
	elseif hasproperty(simulator, :state)
		return simulator.state
	else
		error("Cannot extract state from simulator of type $(typeof(simulator))")
	end
end

"""
	calculate_current_state(sim::BranchedSimulator, path_idx::Int) -> Any

Calculate the current state for a specific path by applying 
all stored instructions to a fresh initial state using the evolve! function.
"""
function calculate_current_state(sim::BranchedSimulator, path_idx::Int)
	# Create a fresh simulator of the same type as the base simulator
	base_sim = similar(sim.base_simulator)
	
	# Set the qubit count to match the current path
	n_qubits = sim.n_qubits[path_idx]
	
	# Initialize the simulator with the right number of qubits
	reinit!(base_sim, n_qubits, 0)
	
	# Apply all instructions using the evolve! function
	evolve!(base_sim, sim.instruction_sequences[path_idx])
	
	# Return the resulting state
	return get_state(base_sim)
end

"""
	pad_bit(i::Int, pos::Int) -> Int

Insert a 0 bit at position pos in the binary representation of i.
This is used for bit manipulation when applying gates to qubits.
"""
function pad_bit(i::Int, pos::Int)
	# Create a mask for bits before the position
	mask_before = (1 << pos) - 1
	# Create a mask for bits after the position
	mask_after = ~mask_before

	# Extract bits before and after the position
	bits_before = i & mask_before
	bits_after = i & mask_after

	# Shift bits after the position to make room for the new bit
	bits_after = bits_after << 1

	# Combine the bits
	return bits_before | bits_after
end

"""
	flip_bit(i::Int, pos::Int) -> Int

Flip the bit at position pos in the binary representation of i.
This is used for bit manipulation when applying gates to qubits.
"""
function flip_bit(i::Int, pos::Int)
	# Create a mask with a 1 at the position to flip
	mask = 1 << pos

	# XOR with the mask to flip the bit
	return i ⊻ mask
end

"""
	allocate_shots(sim::BranchedSimulator) -> Dict{Int, Int}

Allocate shots to each active path based on the number of shots assigned to each path.
Returns a dictionary mapping path indices to shot counts.
"""
function allocate_shots(sim::BranchedSimulator)
	# Create dictionary mapping path indices to shot counts
	return Dict(path_idx => sim.shots[i] for (i, path_idx) in enumerate(sim.active_paths))
end
