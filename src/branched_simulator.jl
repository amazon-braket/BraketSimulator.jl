using Quasar: ClassicalVariable, FunctionDefinition, AbstractGateDefinition

"""
    FramedVariable

A variable type that includes a frame_number field to track the scope frame where the variable was declared.
Similar to ClassicalVariable but with frame tracking.
"""
mutable struct FramedVariable
    # Variable name
    name::String
    
    # Variable type
    type::Any
    
    # Variable value
    val::Any
    
    # Whether the variable is constant
    is_const::Bool
    
    # The frame number where this variable was declared
    # 0 for global scope, 1 for first function call frame, etc.
    frame_number::Int
    
    # Constructor
    function FramedVariable(name::String, type::Any, val::Any, is_const::Bool, frame_number::Int)
        new(name, type, val, is_const, frame_number)
    end
    
    # Constructor from ClassicalVariable
    function FramedVariable(var::ClassicalVariable, frame_number::Int)
        new(var.name, var.type, var.val, var.is_const, frame_number)
    end
end

"""
	BranchedSimulator

A simulator that supports multiple execution paths resulting from mid circuit measurements.
Every time a measurement occurs, all the paths in the active_paths are measured and two new
paths are created depending on whether the measurement was a 0 or 1.

This implementation stores only instruction sequences, not statevectors in order to minimize
the memory requirements.

When a measurement is needed, the state is calculated from scratch by applying
all instructions from the beginning. This state is then used to calculate the probabilities
of measurement for each path. The number of shots passed to each state is then calculated
and any state with 0 shots is then removed from the simulation.

For control sequences that limit the number of paths that enter a set of instructions, the
active_paths list is updated, containing only the indices of the paths that will be evolved
by the instructions
"""
mutable struct BranchedSimulator <: AbstractSimulator
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

	# To denote the current frame that the simulation is working in
	# For the global frame, we will let the frame value be 0.
	curr_frame::Int

	# Classical variable values for each path
	variables::Vector{Dict{String, FramedVariable}}

	# Maps qubit names to their indices in the state vector
	# Maps qubit_name => qubit_index or qubit_name => Vector{Int} for registers
	qubit_mapping::Dict{String, Union{Int, Vector{Int}}}

	# Total number of defined qubits present
	n_qubits::Int

	# Number of shots for each simulation, must be > 0
	shots::Vector{Int}

	# Maps function names to FunctionDefinition objects
	function_defs::Dict{String, FunctionDefinition}

	# Maps built in function names to corresponding Julia function
	function_builtin::Dict{String, Function}

	# Maps gate names to GateDefinition objects
	gate_defs::Dict{String, AbstractGateDefinition}

	# Maps OpenQASM gate names to Julia gate types
	gate_name_mapping::Dict{String, Symbol}

	# Input parameters for the circuit
	inputs::Dict{String, <:Any}

	# Return values for function calls
	return_values::Dict{Int, Any}

	# Simulation indices for continue in for loop
	continue_paths::Vector{Int}

	# Constructor
	function BranchedSimulator(simulator::AbstractSimulator; inputs::Dict{String, <:Any} = Dict{String, Any}())
		# If the number of shots for the simulator is <= 0, then throw an error
		if simulator.shots <= 0
			error("The number of shots for simulating a circuit with MCM must be > 0")
		end

		# Initialize function and gate definitions
		function_defs = Dict{String, FunctionDefinition}()
		gate_defs = Dict{String, AbstractGateDefinition}()

		# Create gate name mapping
		gate_name_mapping = Dict{String, Symbol}(
			# Standard gates with first letter capitalized
			"x" => :X, "y" => :Y, "z" => :Z, "h" => :H,
			"s" => :S, "si" => :Si, "t" => :T, "ti" => :Ti,
			"v" => :V, "vi" => :Vi, "i" => :I,
			"rx" => :Rx, "ry" => :Ry, "rz" => :Rz,
			"cnot" => :CNot, "cy" => :CY, "cz" => :CZ, "cv" => :CV,
			"swap" => :Swap, "iswap" => :ISwap, "pswap" => :PSwap,
			"cswap" => :CSwap, "ccnot" => :CCNot,
			"xx" => :XX, "xy" => :XY, "yy" => :YY, "zz" => :ZZ,
			"ecr" => :ECR, "ms" => :MS,
			"gpi" => :GPi, "gpi2" => :GPi2, "prx" => :PRx,
			"U" => :U,
			# Special cases
			"phaseshift" => :PhaseShift,
			"cphaseshift" => :CPhaseShift,
			"cphaseshift00" => :CPhaseShift00,
			"cphaseshift01" => :CPhaseShift01,
			"cphaseshift10" => :CPhaseShift10,
			"gphase" => :GPhase
		)

		functions_builtin = Dict{String, Function}(
			"arccos" => acos, "arcsin" => asin, "arctan" => atan, "ceiling" => ceil, "cos" => cos, "exp" => exp,
			"floor" => floor, "log" => log, "mod" => mod, "sin" => sin, "sqrt" => sqrt, "tan" => tan
		)

		# Initialize with empty instruction sequence
		instructions = Vector{Instruction}()

		new(
			simulator,
			[instructions],
			[1],
			[Dict{String, Vector{Int}}()],  # Initialize measurements with qubit names as keys
			0,                              # Start with global frame (0)
			[Dict{String, FramedVariable}()], # Initialize variables with framed variables
			Dict{String, Union{Int, Vector{Int}}}(),  # Initialize qubit mapping with qubit names as keys
			0,
			[simulator.shots],
			function_defs,
			functions_builtin,
			gate_defs,
			gate_name_mapping,
			inputs,
			Dict{Int, Any}(),
			Vector{Int}()
		)
	end
end

# Required interface methods
qubit_count(sim::BranchedSimulator) = maximum(sim.n_qubits)
device_id(sim::BranchedSimulator) = device_id(sim.base_simulator)

function reinit!(sim::BranchedSimulator, qubit_count::Int, shots::Int)
	reinit!(sim.base_simulator, qubit_count, shots)
	sim.n_qubits = qubit_count
	sim.shots = [shots]
	sim.active_paths = [1]
	sim.instruction_sequences = [Vector{Instruction}()]
	sim.variables = [Dict{String, FramedVariable}()]
	sim.curr_frame = 0
	sim.gate_defs = Dict{String, AbstractGateDefinition}()
	sim.function_defs = Dict{String, FunctionDefinition}()
	sim.qubit_mapping = Dict{String, Union{Int, Vector{Int}}}()
	sim.measurements = [Dict{String, Vector{Int}}()]
	sim.function_builtin = Dict{String, Function}(
			"arccos" => acos, "arcsin" => asin, "arctan" => atan, "ceiling" => ceil, "cos" => cos, "exp" => exp,
			"floor" => floor, "log" => log, "mod" => mod, "sin" => sin, "sqrt" => sqrt, "tan" => tan
		)
	sim.gate_name_mapping = Dict{String, Symbol}(
			# Standard gates with first letter capitalized
			"x" => :X, "y" => :Y, "z" => :Z, "h" => :H,
			"s" => :S, "si" => :Si, "t" => :T, "ti" => :Ti,
			"v" => :V, "vi" => :Vi, "i" => :I,
			"rx" => :Rx, "ry" => :Ry, "rz" => :Rz,
			"cnot" => :CNot, "cy" => :CY, "cz" => :CZ, "cv" => :CV,
			"swap" => :Swap, "iswap" => :ISwap, "pswap" => :PSwap,
			"cswap" => :CSwap, "ccnot" => :CCNot,
			"xx" => :XX, "xy" => :XY, "yy" => :YY, "zz" => :ZZ,
			"ecr" => :ECR, "ms" => :MS,
			"gpi" => :GPi, "gpi2" => :GPi2, "prx" => :PRx,
			"U" => :U,
			# Special cases
			"phaseshift" => :PhaseShift,
			"cphaseshift" => :CPhaseShift,
			"cphaseshift00" => :CPhaseShift00,
			"cphaseshift01" => :CPhaseShift01,
			"cphaseshift10" => :CPhaseShift10,
			"gphase" => :GPhase
		)

	sim.inputs = Dict{String, Any}()
	sim.return_values = Dict{Int, Any}()
	sim.continue_paths = Vector{Int}()
end


########################
# Measurement Handlers #
########################

"""
	measure_qubit(sim::BranchedSimulator, path_idx::Int, qubit_idx::Int, qubit_name::String)

Starting point for all simulation measurements.

Measure a qubit on a specific path, potentially creating new simulation paths.
Calculates the current state from scratch by applying all instructions.
Returns true if the measurement had two possible outcomes and false otherwise.
The state size is preserved, and half of the states have their probability set to 0.
"""
function measure_qubit(sim::BranchedSimulator, path_idx::Int, qubit_idx::Int, qubit_name::String)
	# Calculate current state and probabilities by applying all instructions from scratch
	current_state = calculate_current_state(sim, path_idx)
	probs = get_measurement_probabilities(current_state, qubit_idx)

	path_shots = sim.shots[path_idx]
	shots_for_outcome_1 = round(Int, path_shots * probs[2])
	shots_for_outcome_0 = path_shots - shots_for_outcome_1

	# We branch only if there is a non zero probability of measuring both outcomes
	if shots_for_outcome_1 > 0 && shots_for_outcome_0 > 0
		handle_branching_measurement!(sim, path_idx, qubit_name, shots_for_outcome_0, shots_for_outcome_1)
		return true
	else
		handle_deterministic_measurement!(sim, path_idx, qubit_name, shots_for_outcome_0 > 0 ? 0 : 1)
		return false
	end
end


"""
	store_measurement_outcome!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, outcome::Int)

Store the measurement outcome for a specific qubit in a specific path.
If the qubit already has measurement outcomes, append the new outcome.
"""
function store_measurement_outcome!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, outcome::Int)
	# Check if the qubit already has measurement outcomes
	if haskey(sim.measurements[path_idx], qubit_name)
		push!(sim.measurements[path_idx][qubit_name], outcome)
	else
		sim.measurements[path_idx][qubit_name] = [outcome]
	end
end

"""
	handle_branching_measurement!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, probs::Vector{Float64})

Handle the case where measurement can result in multiple outcomes, creating a new branch.
Keep the full state size and set probabilities to 0 for the states that don't match the measurement outcome.
"""
function handle_branching_measurement!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, shots_for_outcome_0::Int, shots_for_outcome_1::Int)

	sim.shots[path_idx] = shots_for_outcome_0

	# Create new path for outcome 1
	# Append the path at the end of the current states list
	new_path_idx = length(sim.instruction_sequences) + 1
	push!(sim.shots, shots_for_outcome_1)

	# Deep copy measurements to avoid aliasing between paths
	new_measurements = Dict{String, Vector{Int}}()
	for (qubit_name, outcomes) in sim.measurements[path_idx]
		new_measurements[qubit_name] = copy(outcomes)
	end
	push!(sim.measurements, new_measurements)
	new_operations = copy(sim.instruction_sequences[path_idx])
	push!(sim.instruction_sequences, new_operations)

	# Deep copy the framed variables to avoid aliasing between paths
	new_variables = Dict{String, FramedVariable}()
	for (var_name, var) in sim.variables[path_idx]
		new_val = var.val isa Vector ? copy(var.val) : var.val
		new_variables[var_name] = FramedVariable(var.name, var.type, new_val, var.is_const, var.frame_number)
	end
	push!(sim.variables, new_variables)

	# Store the measurement outcome for both paths
	store_measurement_outcome!(sim, new_path_idx, qubit_name, 1)
	store_measurement_outcome!(sim, path_idx, qubit_name, 0)
end

"""
	handle_deterministic_measurement!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, probs::Vector{Float64})

Handle the case where measurement has only one possible outcome. No need for branching or anything, just store
the value of the classical variables.
"""
function handle_deterministic_measurement!(sim::BranchedSimulator, path_idx::Int, qubit_name::String, outcome::Int)
	# Store the measurement outcome
	store_measurement_outcome!(sim, path_idx, qubit_name, outcome)
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
	set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, value)

Set a classical variable value for a specific path.
"""
function set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, value::ClassicalVariable)
	# Convert ClassicalVariable to FramedVariable with current frame number
	sim.variables[path_idx][var_name] = FramedVariable(value, sim.curr_frame)
end

"""
	set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, type, val, is_const)

Set a variable value for a specific path with individual components.
"""
function set_variable!(sim::BranchedSimulator, path_idx::Int, var_name::String, type::Any, val::Any, is_const::Bool)
	sim.variables[path_idx][var_name] = FramedVariable(var_name, type, val, is_const, sim.curr_frame)
end

"""
	get_variable(sim::BranchedSimulator, path_idx::Int, var_name::String, default=nothing)

Get a classical variable value for a specific path.
"""
function get_variable(sim::BranchedSimulator, path_idx::Int, var_name::String, default = nothing)
	return haskey(sim.variables[path_idx], var_name) ? sim.variables[path_idx][var_name] : default
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
	n_qubits = sim.n_qubits
	
	# Make sure n_qubits is at least 1
	if n_qubits <= 0
		n_qubits = 1
	end

	# Initialize the simulator with the right number of qubits and shots
	reinit!(base_sim, n_qubits, sim.shots[path_idx])

	# Apply all instructions using the evolve! function
	# TODO: Add parallelization here
	evolve!(base_sim, sim.instruction_sequences[path_idx])

	# Debug: Print the final state
	final_state = get_state(base_sim)

	# Return the resulting state
	return final_state
end

"""
	calculate_current_state(sim::BranchedSimulator) -> Dict{Int, Any}

Calculate the current state for all active paths by applying 
all stored instructions to fresh initial states using the evolve! function.
Returns a dictionary mapping path indices to their respective states.
"""
function calculate_current_state(sim::BranchedSimulator)
	# Create a dictionary to store results for each path
	results = Dict{Int, Any}()
	
	# Calculate state for each active path
	for path_idx in sim.active_paths
		path_state = calculate_current_state(sim, path_idx)
		results[path_idx] = path_state
	end
	
	return results
end
