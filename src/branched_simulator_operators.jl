using Quasar: ClassicalVariable, FunctionDefinition, AbstractGateDefinition

"""
	BranchedSimulatorOperators

A simulator that supports multiple execution paths resulting from measurements.
This implementation stores only instruction sequences, not statevectors.
When a measurement is needed, the state is calculated from scratch by applying
all instructions from the beginning.
"""
mutable struct BranchedSimulatorOperators <: AbstractSimulator
	# Base simulator type
	base_simulator::AbstractSimulator

	# Collection of instruction sequences for different paths
	# Each path stores a list of instructions to apply
	instruction_sequences::Vector{Vector{Instruction{<:Operator}}}
	
	# Active paths (indices into instruction_sequences that are still evolving)
	active_paths::Vector{Int}

	# Tracks and centralizes the measurement outcomes of qubits
	# Maps path_idx => Dict(qubit_name => measurement_outcomes)
	measurements::Vector{Dict{String, Vector{Int}}}

	# Classical variable values for each path
	variables::Vector{Dict{String, ClassicalVariable}}

	# Maps qubit names to their indices in the state vector
	# Maps path_idx => Dict(qubit_name => qubit_index)
	qubit_mapping::Dict{String, Int}

	# Total number of defined qubits present
	n_qubits::Int

	# Number of shots for each simulation
	shots::Vector{Int}

	# Maps function names to FunctionDefinition objects
	function_defs::Dict{String, FunctionDefinition}

	# Maps gate names to GateDefinition objects
	gate_defs::Dict{String, AbstractGateDefinition}

	# Input parameters for the circuit
	inputs::Dict{String, Any}

	# Constructor
	function BranchedSimulatorOperators(simulator::AbstractSimulator; inputs::Dict{String, Any} = Dict{String, Any}())
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
		
		# Initialize with empty instruction sequence
		instructions = Vector{Instruction}()

		new(
			simulator,
			[instructions],  # Initialize with empty instruction sequence
			[1],
			[Dict{String, Vector{Int}}()],  # Initialize measurements with qubit names as keys
			[Dict{String, ClassicalVariable}()],
			Dict{String, Int}(),          # Initialize qubit mapping with qubit names as keys
			0,
			[simulator.shots],
			function_defs,
			gate_defs,
			inputs,
		)
	end
end

# Required interface methods
qubit_count(sim::BranchedSimulatorOperators) = maximum(sim.n_qubits)
device_id(sim::BranchedSimulatorOperators) = device_id(sim.base_simulator)

"""
	store_measurement_outcome!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String, outcome::Int)

Store the measurement outcome for a specific qubit in a specific path.
If the qubit already has measurement outcomes, append the new outcome.
"""
function store_measurement_outcome!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String, outcome::Int)
	# Check if the qubit already has measurement outcomes
	if haskey(sim.measurements[path_idx], qubit_name)
		# Append the new outcome
		push!(sim.measurements[path_idx][qubit_name], outcome)
	else
		# Create a new vector with the outcome
		sim.measurements[path_idx][qubit_name] = [outcome]
	end
end

"""
	handle_branching_measurement!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, 
								 qubit_name::String, probs::Vector{Float64})

Handle the case where measurement can result in multiple outcomes, creating a new branch.
Keep the full state size and set probabilities to 0 for the states that don't match the measurement outcome.
"""
function handle_branching_measurement!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int,
	qubit_name::String, probs::Vector{Float64})

	# Calculate shots for each path based on probabilities
	path_shots = sim.shots[path_idx]
	shots_for_outcome_1 = round(Int, path_shots * probs[2])
	shots_for_outcome_0 = path_shots - shots_for_outcome_1

	# Create new path for outcome 1
	new_path_idx = length(sim.instruction_sequences) + 1
	push!(sim.shots, shots_for_outcome_1)

	# Copy measurements, operations, and variables for the new path
	new_measurements = copy(sim.measurements[path_idx])
	push!(sim.measurements, new_measurements)
	new_operations = copy(sim.instruction_sequences[path_idx])
	push!(sim.instruction_sequences, new_operations)
	
	# Deep copy the variables to avoid aliasing between paths
	new_variables = Dict{String, ClassicalVariable}()
	for (var_name, var) in sim.variables[path_idx]
		# Create a new ClassicalVariable with the same properties but a copy of the value
		new_val = var.val isa Vector ? copy(var.val) : var.val
		new_variables[var_name] = ClassicalVariable(var.name, var.type, new_val, var.is_const)
	end
	push!(sim.variables, new_variables)

	# Store the measurement outcome for the new path
	store_measurement_outcome!(sim, new_path_idx, qubit_name, 1)

	# Update original path for outcome 0
	sim.shots[path_idx] = shots_for_outcome_0

	# Store the measurement outcome for the original path
	store_measurement_outcome!(sim, path_idx, qubit_name, 0)
end

"""
	handle_deterministic_measurement!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, 
									 qubit_name::String, probs::Vector{Float64})

Handle the case where measurement has only one possible outcome.
Keep the full state size and set probabilities to 0 for the states that don't match the measurement outcome.
"""
function handle_deterministic_measurement!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int,
	qubit_name::String, probs::Vector{Float64})
	# Determine the only possible outcome
	outcome = probs[2] > 0.5 ? 1 : 0

	# Store the measurement outcome
	store_measurement_outcome!(sim, path_idx, qubit_name, outcome)
end

"""
	measure_qubit(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, qubit_name::String)

Measure a qubit on a specific path, potentially creating new simulation paths.
Calculates the current state from scratch by applying all instructions.
Returns true if the measurement had two possible outcomes and false otherwise.
The state size is preserved, and half of the states have their probability set to 0.
"""
function measure_qubit(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, qubit_name::String)
	# Calculate current state by applying all instructions from scratch
	current_state = calculate_current_state(sim, path_idx)
	
	# Get measurement probabilities for this path
	probs = get_measurement_probabilities(current_state, qubit_idx)

	# Determine if we need to branch (both outcomes possible) or not
	if probs[1] > 0 && probs[2] > 0
		# Both outcomes are possible - create a branch
		handle_branching_measurement!(sim, path_idx, qubit_idx, qubit_name, probs)
		return true
	else
		# Only one outcome is possible - no branching needed
		handle_deterministic_measurement!(sim, path_idx, qubit_idx, qubit_name, probs)
		return false
	end
end

"""
	measure_qubit!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, qubit_name::String)

Measure a qubit on a specific path, creating a temporary new_active_paths vector.
This is a convenience method that creates a new vector for new_active_paths.
The state size is preserved, and half of the states have their probability set to 0.
"""
function measure_qubit!(sim::BranchedSimulatorOperators, path_idx::Int, qubit_idx::Int, qubit_name::String)
	# We don't need to check if the qubit was previously measured
	# since we allow multiple measurements of the same qubit
	
	# Measure the qubit
	path_split = measure_qubit(sim, path_idx, qubit_idx, qubit_name)
	
	# If this is called directly, we need to update the active paths
	if path_split && path_idx in sim.active_paths
		# Remove the original path
		filter!(p -> p != path_idx, sim.active_paths)
		# Add the new paths
		push!(sim.active_paths, path_idx)
		push!(sim.active_paths, length(sim.instruction_sequences))
	end
	
	return nothing
end

"""
	filter_paths!(sim::BranchedSimulatorOperators, condition_func)

Filter active paths based on a condition function.
"""
function filter_paths!(sim::BranchedSimulatorOperators, condition_func)
	sim.active_paths = filter(path_idx -> condition_func(path_idx, sim), sim.active_paths)
	return nothing
end

"""
	get_measurement_result(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String)

Get the measurement results for a specific qubit in a specific path.
Returns the most recent measurement outcome, or -1 if the qubit has not been measured.
"""
function get_measurement_result(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String)
	# Check if the qubit has been measured
	if haskey(sim.measurements[path_idx], qubit_name)
		# Return the most recent measurement outcome
		return sim.measurements[path_idx][qubit_name][end]
	else
		# Qubit has not been measured
		return -1
	end
end

"""
	get_measurement_results(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String)

Get all measurement results for a specific qubit in a specific path.
Returns an empty vector if the qubit has not been measured.
"""
function get_measurement_results(sim::BranchedSimulatorOperators, path_idx::Int, qubit_name::String)
	# Return all measurement outcomes
	return get(sim.measurements[path_idx], qubit_name, Int[])
end

"""
	set_variable!(sim::BranchedSimulatorOperators, path_idx::Int, var_name::String, value)

Set a classical variable value for a specific path.
"""
function set_variable!(sim::BranchedSimulatorOperators, path_idx::Int, var_name::String, value::ClassicalVariable)
	sim.variables[path_idx][var_name] = value
end

"""
	get_variable(sim::BranchedSimulatorOperators, path_idx::Int, var_name::String, default=nothing)

Get a classical variable value for a specific path.
"""
function get_variable(sim::BranchedSimulatorOperators, path_idx::Int, var_name::String, default = nothing)
	return haskey(sim.variables[path_idx], var_name) ? sim.variables[path_idx][var_name] : default
end


"""
	calculate_current_state(sim::BranchedSimulatorOperators, path_idx::Int) -> Any

Calculate the current state for a specific path by applying 
all stored instructions to a fresh initial state using the evolve! function.
"""
function calculate_current_state(sim::BranchedSimulatorOperators, path_idx::Int)
	# Create a fresh simulator of the same type as the base simulator
	base_sim = similar(sim.base_simulator)
	
	# Set the qubit count to match the current path
	n_qubits = sim.n_qubits
	
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

	flip_bit(i::Int, pos::Int) -> Int

Flip the bit at position pos in the binary representation of i.
This is used for bit manipulation when applying gates to qubits.

function flip_bit(i::Int, pos::Int)
	# Create a mask with a 1 at the position to flip
	mask = 1 << pos

	# XOR with the mask to flip the bit
	return i ⊻ mask
end


	allocate_shots(sim::BranchedSimulatorOperators) -> Dict{Int, Int}

Allocate shots to each active path based on the number of shots assigned to each path.
Returns a dictionary mapping path indices to shot counts.

function allocate_shots(sim::BranchedSimulatorOperators)
	# Create dictionary mapping path indices to shot counts
	return Dict(path_idx => sim.shots[i] for (i, path_idx) in enumerate(sim.active_paths))
end
"""
