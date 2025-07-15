module BraketSimulator

using
	Dates,
	Combinatorics,
	LinearAlgebra,
	JSON3,
	StaticArrays,
	StatsBase,
	UUIDs,
	StructTypes,
	Random,
	Quasar,
	PrecompileTools

export StateVectorSimulator, DensityMatrixSimulator, BranchedSimulator, evolve!, simulate, ValidationError,
	calculate_probability_threshold, prune_paths!, allocate_shots, get_measurement_probabilities, apply_projection!,
	_apply_reset!, measure_qubit!, get_measurement_result, set_variable!, get_variable, qubit_count, device_id, apply_gate!,
	new_to_circuit, _evolve_branched_ast, BranchedSimulatorOperators, evolve_branched_operators, calculate_shots_per_state

const AbstractStateVector{T} = AbstractVector{T}
const AbstractDensityMatrix{T} = AbstractMatrix{T}

abstract type AbstractSimulator end
ap_size(shots::Int, qubit_count::Int) = (shots > 0 && qubit_count < 30) ? 2^qubit_count : 0

"""
	FreeParameter
	FreeParameter(name::Symbol) -> FreeParameter

Struct representing a free parameter, which may be used to initialize
to a parametrized [`Gate`](@ref) or [`Noise`](@ref) and then given a
fixed value later by supplying a mapping to a [`Circuit`](@ref).
"""
struct FreeParameter
	name::Symbol
	FreeParameter(name::Symbol) = new(name)
	FreeParameter(name::String) = new(Symbol(name))
end
Base.copy(fp::FreeParameter) = fp
Base.show(io::IO, fp::FreeParameter) = print(io, string(fp.name))

function complex_matrix_from_ir(mat::Vector{Vector{Vector{T}}}) where {T <: Number}
	m = zeros(complex(T), length(mat), length(mat))
	for ii in 1:length(mat), jj in 1:length(mat)
		m[ii, jj] = complex(mat[ii][jj][1], mat[ii][jj][2])
	end
	return m
end

function complex_matrix_to_ir(m::Matrix{T}) where {T <: Complex}
	mat = Vector{Vector{Vector{real(T)}}}(undef, size(m, 1))
	for row in 1:size(m, 1)
		mat[row] = Vector{Vector{T}}(undef, size(m, 2))
		for col in 1:size(m, 2)
			mat[row][col] = [real(m[row, col]), imag(m[row, col])]
		end
	end
	return mat
end
complex_matrix_to_ir(m) = m


include("raw_schema.jl")
include("qubit_set.jl")
include("operators.jl")
include("observables.jl")
using .Observables
include("results.jl")
include("pragmas.jl")
include("builtin_gates.jl")
include("gates.jl")
include("noises.jl")
include("schemas.jl")
include("circuit.jl")
include("validation.jl")
include("custom_gates.jl")
include("pow_gates.jl")
include("gate_kernels.jl")
include("noise_kernels.jl")
include("result_types.jl")
include("properties.jl")
include("sv_simulator.jl")
include("dm_simulator.jl")
include("branched_simulator.jl")
include("branched_simulator_operators.jl")
include("branched_evolve_operators.jl")
include("branched_evolve.jl")

function __init__()
	Quasar.builtin_gates[] = builtin_gates
	Quasar.parse_pragma[]  = parse_pragma
	Quasar.visit_pragma[]  = visit_pragma
end

const LOG2_CHUNK_SIZE = 16
const CHUNK_SIZE      = 2^LOG2_CHUNK_SIZE

function _index_to_endian_bits(ix::Int, qubit_count::Int)
	bits = Vector{Int}(undef, qubit_count)
	for qubit ∈ 0:(qubit_count-1)
		bit = (ix >> qubit) & 1
		bits[end-qubit] = bit
	end
	return bits
end

function _formatted_measurements(simulator::D, measured_qubits::Vector{Int}) where {D <: AbstractSimulator}
	sim_samples = samples(simulator)
	n_qubits    = qubit_count(simulator)
	return [_index_to_endian_bits(sample, n_qubits)[measured_qubits .+ 1] for sample in sim_samples]
end

function _build_metadata(simulator, ir; shots = nothing)
	task_mtd = TaskMetadata(
		braketSchemaHeader("braket.task_result.task_metadata", "1"),
		string(uuid4()),
		shots === nothing ? simulator.shots : shots,
		device_id(simulator),
		nothing,
		nothing,
		nothing,
		nothing,
		nothing,
	)
	addl_mtd = AdditionalMetadata(
		ir,
		nothing,
		nothing,
		nothing,
		nothing,
		nothing,
		nothing,
		nothing,
	)
	return task_mtd, addl_mtd
end

function _bundle_results(
	results::Vector{ResultTypeValue},
	circuit_ir,
	simulator::D,
	measured_qubits = Set{Int}(0:(qubit_count(simulator)-1)),
) where {D <: AbstractSimulator}
	sorted_qubits = sort(collect(measured_qubits))
	formatted_samples = if simulator.shots > 0
		_formatted_measurements(simulator, sorted_qubits)
	else
		nothing
	end
	return GateModelTaskResult(
		braketSchemaHeader("braket.task_result.gate_model_task_result", "1"),
		formatted_samples,
		nothing,
		results,
		sorted_qubits,
		_build_metadata(simulator, circuit_ir)...,
	)
end

function _generate_results(
	result_types,
	simulator::D,
) where {D <: AbstractSimulator}
	result_values = map(result_type -> calculate(result_type, simulator), result_types)
	ir_results    = map(StructTypes.lower, result_types)
	return [ResultTypeValue(ir_results[r_ix], complex_matrix_to_ir(result_values[r_ix])) for r_ix in 1:length(result_values)]
end

_translate_result_type(r::IR.Amplitude)     = Amplitude(r.states)
_translate_result_type(r::IR.StateVector)   = StateVector()
_translate_result_type(r::IR.DensityMatrix) = DensityMatrix(r.targets)
_translate_result_type(r::IR.Probability)   = Probability(r.targets)
for (RT, IRT) in ((:Expectation, :(IR.Expectation)), (:Variance, :(IR.Variance)), (:Sample, :(IR.Sample)))
	@eval begin
		function _translate_result_type(r::$IRT)
			obs = StructTypes.constructfrom(Observables.Observable, r.observable)
			$RT(obs, QubitSet(r.targets))
		end
	end
end
_translate_result_types(results::Vector{AbstractProgramResult}) = map(_translate_result_type, results)

function _compute_exact_results(d::AbstractSimulator, program::Program, qubit_count::Int)
	result_types = _translate_result_types(program.results)
	_validate_result_types_qubits_exist(result_types, qubit_count)
	return _generate_results(result_types, d)
end

function _compute_exact_results(d::AbstractSimulator, program::Circuit, qubit_count::Int)
	_validate_result_types_qubits_exist(program.result_types, qubit_count)
	return _generate_results(program.result_types, d)
end

"""
	_get_measured_qubits(program, qubit_count::Int) -> Vector{Int}

Get the qubits measured by the program. If [`Measure`](@ref)
instructions are present in the program's instruction list,
their targets are used to form the list of measured qubits.
If not, all qubits from 0 to `qubit_count-1` are measured. 
"""
function _get_measured_qubits(program, qubit_count::Int)
	measure_inds = findall(ix->ix.operator isa Measure, program.instructions)
	isempty(measure_inds) && return collect(0:(qubit_count-1))
	measure_ixs     = program.instructions[measure_inds]
	measure_targets = (convert(Vector{Int}, measure.target) for measure in measure_ixs)
	measured_qubits = unique(reduce(vcat, measure_targets, init = Int[]))
	return measured_qubits
end
"""
	_prepare_program(circuit_ir::OpenQasmProgram, inputs::Dict{String, <:Any}, shots::Int) -> (Program, Int)

Parse the OpenQASM3 source, apply any `inputs` provided for the simulation, and compute
basis rotation instructions if running with non-zero shots. Return the `Program` after
parsing and the qubit count of the circuit.
"""
function _prepare_program(circuit_ir::OpenQasmProgram, inputs::Dict{String, <:Any}, shots::Int)
	ir_inputs     = isnothing(circuit_ir.inputs) ? Dict{String, Float64}() : circuit_ir.inputs
	merged_inputs = merge(ir_inputs, inputs)
	src           = circuit_ir.source::String
	circuit       = to_circuit(src, merged_inputs)
	n_qubits      = qubit_count(circuit)
	if shots > 0
		_verify_openqasm_shots_observables(circuit, n_qubits)
		basis_rotation_instructions!(circuit)
	end
	return circuit, n_qubits
end

"""
	_combine_operations(program, shots::Int) -> Program

Combine explicit instructions and basis rotation instructions (if necessary).
Validate that all operations are performed on qubits within `qubit_count`.
"""
function _combine_operations(program, shots::Int)
	operations = program.instructions
	if shots > 0 && !isempty(program.basis_rotation_instructions)
		operations = vcat(operations, program.basis_rotation_instructions)
	end
	_validate_operation_qubits(operations)
	return operations
end
"""
	_compute_results(::, simulator, program::Program, n_qubits::Int, shots::Int) -> Vector{ResultTypeValue}

Compute the results once `simulator` has finished applying all the instructions. The results depend on the IR type if `shots>0`:

- For JAQCD IR (`Program`), the results array is *empty* because the Braket SDK computes the results from the IR directly.
- For OpenQASM IR (`OpenQasmProgram`), the results array is *empty* only if no results are present in the parsed IR. Otherwise,
  the results array is populated with the parsed result types (to help the Braket SDK compute them from the sampled measurements)
  and a placeholder zero value.
"""
function _compute_results(simulator, program::Circuit, n_qubits, shots)
	results          = program.result_types
	has_no_results   = isnothing(results) || isempty(results)
	analytic_results = shots == 0 && !has_no_results
	if analytic_results
		return _compute_exact_results(simulator, program, n_qubits)
	elseif has_no_results
		return ResultTypeValue[]
	else
		return ResultTypeValue[ResultTypeValue(StructTypes.lower(result_type), 0.0) for result_type in results]
	end
end

function _validate_circuit_ir(simulator, circuit_ir::Circuit, qubit_count::Int, shots::Int)
	_validate_ir_results_compatibility(simulator, circuit_ir.result_types, Val(:JAQCD))
	_validate_ir_instructions_compatibility(simulator, circuit_ir, Val(:JAQCD))
	_validate_shots_and_ir_results(shots, circuit_ir.result_types, qubit_count)
	return
end

"""
	simulate(simulator::AbstractSimulator, circuit_ir::Union{OpenQasmProgram, Program}, shots::Int; kwargs...) -> GateModelTaskResult

Simulate the evolution of a state vector or density matrix using the passed-in `simulator`.
The instructions to apply (gates and noise channels) and measurements to make are
encoded in `circuit_ir`. Supported IR formats are `OpenQASMProgram` (OpenQASM3)
and `Program` (JAQCD). Returns a `GateModelTaskResult` containing the individual shot
measurements (if `shots > 0`), final calculated results, circuit IR, and metadata
about the task.
"""
function simulate(
	simulator::AbstractSimulator,
	circuit_ir::OpenQasmProgram,
	shots::Int;
	inputs = Dict{String, Float64}(),
	kwargs...,
)
	program, n_qubits = _prepare_program(circuit_ir, inputs, shots)
	_validate_circuit_ir(simulator, program, n_qubits, shots)
	operations = _combine_operations(program, shots)
	reinit!(simulator, n_qubits, shots)
	simulator       = evolve!(simulator, operations)
	measured_qubits = _get_measured_qubits(program, n_qubits)
	results         = _compute_results(simulator, program, n_qubits, shots)
	return _bundle_results(results, circuit_ir, simulator, measured_qubits)
end


"""
	_prepare_program(circuit_ir::OpenQasmProgram, inputs::Dict{String, <:Any}, shots::Int) -> (QasmExpression, Int)

Parse the OpenQASM3 source and compute basis rotation instructions if running with non-zero shots. 
Return the `Program` after parsing and the qubit count of the circuit.
"""
function _new_prepare_program(circuit_ir::OpenQasmProgram, shots::Int)
	src      = circuit_ir.source::String
	circuit  = new_to_circuit(src)
	n_qubits = qubit_count(circuit)
	if shots > 0
		_verify_openqasm_shots_observables(circuit, n_qubits)
		basis_rotation_instructions!(circuit)
	end
	return circuit, n_qubits
end

# new_to_circuit is already defined in circuit.jl

"""
	simulate(simulator::AbstractSimulator, circuit_irs::Vector{<:Union{Program, OpenQasmProgram}}, shots::Int; max_parallel::Int=min(32, Threads.nthreads()), inputs=Dict{String,Float64}(), kwargs...) -> Vector{GateModelTaskResult}

Simulate the evolution of a *batch* of state vectors or density matrices using the passed in `simulator`.
The instructions to apply (gates and noise channels) and measurements to make are
encoded in `circuit_irs`. Supported IR formats are `OpenQASMProgram` (OpenQASM3)
and `Program` (JAQCD).

The simulation of the batch is done in parallel using threads.
The keyword argument `max_parallel` specifies the number of evolutions to simulate in
parallel -- the default value is whichever of `32` and `Threads.nthreads()` is smaller.
This is to avoid overwhelming the thread scheduler with too many small tasks waiting to
run, as each evolution is itself threaded. This value may change in the future.

The `inputs` keyword argument can be a `Dict{String}` or a `Vector{Dict{String}}`. In
the first case, the same input values are applied to all `circuit_irs`. In the second,
the length of the `inputs` *must* be the same as the length of `circuit_irs`, and the
`n`-th `inputs` is applied to the `n`-th `circuit_irs`.

Returns a `Vector{GateModelTaskResult}`, each element of which contains the individual shot
measurements (if `shots > 0`), final calculated results, corresponding circuit IR, and metadata
about the task.
"""
function simulate(simulator::AbstractSimulator,
	circuit_irs::Vector{<:Union{Program, OpenQasmProgram}},
	shots::Int = 0;
	max_parallel::Int = min(32, Threads.nthreads()),
	inputs = Dict{String, Float64}(),
	kwargs...,
)

	n_tasks         = length(circuit_irs)
	is_single_task  = n_tasks == 1
	is_single_input = inputs isa Dict || length(inputs) == 1
	if is_single_input
		single_input = inputs isa Vector ? only(inputs) : inputs
		if is_single_task
			return [simulate(simulator, only(circuit_irs), shots; inputs = single_input, kwargs...)]
		else
			inputs = [deepcopy(single_input) for ix in 1:n_tasks]
		end
	end
	!is_single_task && !is_single_input && (n_tasks != length(inputs)) && throw(ArgumentError("number of inputs ($(length(inputs))), and circuit IRs ($n_tasks) must be equal."))
	todo_tasks_ch = Channel{Int}(ch->foreach(ix->put!(ch, ix), 1:n_tasks), n_tasks)
	n_task_threads = min(max_parallel, n_tasks)
	results = Vector{GateModelTaskResult}(undef, n_tasks)
	function process_work(my_sim)
		while isready(todo_tasks_ch)
			my_ix = -1
			# need to lock the channel as it may become empty
			# and "unready" in between the while-loop call
			# and the call to take!
			lock(todo_tasks_ch) do
				my_ix = isready(todo_tasks_ch) ? take!(todo_tasks_ch) : -1
			end
			# if my_ix is still -1, the channel is empty and
			# there's no more work to do
			my_ix == -1 && break
			results[my_ix] = simulate(my_sim, circuit_irs[my_ix], shots; inputs = inputs[my_ix])
		end
		return
	end
	# need sync here to ensure all the spawns launch
	tasks = @sync [Threads.@spawn process_work(similar(simulator)) for worker in 1:n_task_threads]
	# tasks don't return anything so we can wait rather than fetch
	wait.(tasks)
	# check to ensure all the results were in fact populated
	for r_ix in 1:n_tasks
		@assert isassigned(results, r_ix)
	end
	return results
end


function _mmap_large_result_values(results)
	to_mmap = findall(rt->sizeof(rt.value) > 2^20, results.resultTypes)
	isempty(to_mmap) && return nothing, nothing
	mmap_files  = String[]
	obj_lengths = Int[]
	for r_ix in to_mmap
		push!(obj_lengths, length(results.resultTypes[r_ix].value))
		tmp_path, io = mktemp()
		write(io, results.resultTypes[r_ix].value)
		empty!(results.resultTypes[r_ix].value)
		close(io)
		push!(mmap_files, tmp_path)
	end
	py_paths    = tuple(mmap_files...)
	py_lens     = tuple(obj_lengths...)
	mmap_files  = nothing
	obj_lengths = nothing
	return py_paths, py_lens
end

function BraketSimulator.simulate(simulator_id::String, task_spec::String, py_inputs::String, shots::Int; kwargs...)
	inputs = JSON3.read(py_inputs, Dict{String, Any})
	jl_spec = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, inputs)
	simulator = create_sim(simulator_id, shots)
	jl_results = simulate(simulator, jl_spec, shots; kwargs...)
	py_paths, py_lens = _mmap_large_result_values(jl_results)
	py_results = JSON3.write(jl_results)
	simulator = nothing
	inputs = nothing
	jl_spec = nothing
	jl_results = nothing
	return py_results, py_paths, py_lens
end
function BraketSimulator.simulate(simulator_id::String, task_specs::AbstractVector, py_inputs::String, shots::Int; kwargs...)
	inputs = JSON3.read(py_inputs, Vector{Dict{String, Any}})
	jl_specs = map(zip(task_specs, inputs)) do (task_spec, input)
		BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, input)
	end
	simulator = create_sim(simulator_id, shots)
	jl_results = simulate(simulator, jl_specs, shots; kwargs...)
	paths_and_lens = JSON3.write(map(_mmap_large_result_values, jl_results))
	jsons = JSON3.write(jl_results)
	simulator = nothing
	jl_results = nothing
	inputs = nothing
	jl_specs = nothing
	return jsons, paths_and_lens
end


"""
	simulate(simulator::AbstractSimulator, circuit_ir::Union{OpenQasmProgram, Program}, shots::Int; kwargs...) -> GateModelTaskResult

Simulate the evolution of a state vector or density matrix using the passed-in `simulator`.
The instructions to apply (gates and noise channels) and measurements to make are
encoded in `circuit_ir`. Supported IR format is `OpenQASMProgram`. Returns a list of measured qubit
results, combining both MCM and terminal measurements.
"""

function simulate(
	simulator::BranchedSimulatorOperators,
	circuit_ir::OpenQasmProgram,
	shots::Int;
	inputs = Dict{String, Float64}(),
	kwargs...,
)
	if shots <= 0
		error("The number of inputted shots must be a positive integer")
	end
	if shots != simulator.shots[1]
		error("The number of shots in the simulator must be equal to the number of shots passed in")
	end

	# Parse and prepare the program
	program = new_to_circuit(circuit_ir.source)
	# _validate_circuit_ir(simulator, program, n_qubits, shots)

	# Initialize the simulator
	reinit!(simulator, 0, shots) # Fix this function

	# Use the branched simulator to handle measurements and control flow
	branched_sim = evolve_branched_operators(simulator, program, inputs)

	circuit_measurements = calculate_shots_per_state(branched_sim)
	# Create the result object
	measured_qubits = 0:(branched_sim.n_qubits-1)

	# Return the result
	return GateModelTaskResult(
		braketSchemaHeader("braket.task_result.gate_model_task_result", "1"),
		circuit_measurements,
		nothing,
		ResultTypeValue[],  # For now, no result types
		measured_qubits,
		_build_metadata(branched_sim, circuit_ir; shots = branched_sim.shots[1])...,
	)
end

# these functions are for calls from an "external" language
# like Python or Rust, when we're calling this package from a
# separate process and thus don't want to have to IPC large
# blobs of data back and forth/deal with having to serialize
# Julia objects
function create_sim(simulator_id::String, shots::Int)
	return if simulator_id == "braket_sv_v2"
		StateVectorSimulator(0, shots)
	elseif simulator_id == "braket_dm_v2"
		DensityMatrixSimulator(0, shots)
	elseif simulator_id == "braket_sv_branched"
		BranchedSimulatorOperators(StateVectorSimulator(0, shots)) # Fix this
	end
end

@setup_workload begin
	custom_qasm = """
			   int[8] two = 2;
			   gate x a { U(π, 0, π) a; }
			   gate cx c, a {
				   pow(1) @ ctrl @ x c, a;
			   }
			   gate cxx_1 c, a {
				   pow(two) @ cx c, a;
			   }
			   gate cxx_2 c, a {
				   pow(1/2) @ pow(4) @ cx c, a;
			   }
			   gate cxxx c, a {
				   pow(1) @ pow(two) @ cx c, a;
			   }

			   qubit q1;
			   qubit q2;
			   qubit q3;
			   qubit q4;
			   qubit q5;

			   pow(1/2) @ x q1;       // half flip
			   pow(1/2) @ x q1;       // half flip
			   cx q1, q2;   // flip
			   cxx_1 q1, q3;    // don't flip
			   cxx_2 q1, q4;    // don't flip
			   cnot q1, q5;    // flip
			   x q3;       // flip
			   x q4;       // flip

			   s q1;   // sqrt z
			   s q1;   // again
			   inv @ z q1; // inv z
			   """;
	noise_qasm = """
		qubit[2] qs;

		#pragma braket noise bit_flip(.5) qs[1]
		#pragma braket noise phase_flip(.5) qs[0]
		#pragma braket noise pauli_channel(.1, .2, .3) qs[0]
		#pragma braket noise depolarizing(.5) qs[0]
		#pragma braket noise two_qubit_depolarizing(.9) qs
		#pragma braket noise two_qubit_depolarizing(.7) qs[1], qs[0]
		#pragma braket noise two_qubit_dephasing(.6) qs
		#pragma braket noise amplitude_damping(.2) qs[0]
		#pragma braket noise generalized_amplitude_damping(.2, .3)  qs[1]
		#pragma braket noise phase_damping(.4) qs[0]
		#pragma braket noise kraus([[0.9486833im, 0], [0, 0.9486833im]], [[0, 0.31622777], [0.31622777, 0]]) qs[0]
		#pragma braket noise kraus([[0.9486832980505138, 0, 0, 0], [0, 0.9486832980505138, 0, 0], [0, 0, 0.9486832980505138, 0], [0, 0, 0, 0.9486832980505138]], [[0, 0.31622776601683794, 0, 0], [0.31622776601683794, 0, 0, 0], [0, 0, 0, 0.31622776601683794], [0, 0, 0.31622776601683794, 0]]) qs[{1, 0}]
		"""
	unitary_qasm = """
		qubit[3] q;

		x q[0];
		h q[1];

		// unitary pragma for t gate
		#pragma braket unitary([[1.0, 0], [0, 0.70710678 + 0.70710678im]]) q[0]
		ti q[0];

		// unitary pragma for h gate (with phase shift)
		#pragma braket unitary([[0.70710678im, 0.70710678im], [0 - -0.70710678im, -0.0 - 0.70710678im]]) q[1]
		gphase(-π/2) q[1];
		h q[1];

		// unitary pragma for ccnot gate
		#pragma braket unitary([[1.0, 0, 0, 0, 0, 0, 0, 0], [0, 1.0, 0, 0, 0, 0, 0, 0], [0, 0, 1.0, 0, 0, 0, 0, 0], [0, 0, 0, 1.0, 0, 0, 0, 0], [0, 0, 0, 0, 1.0, 0, 0, 0], [0, 0, 0, 0, 0, 1.0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1.0], [0, 0, 0, 0, 0, 0, 1.0, 0]]) q
		"""
	sv_adder_qasm = """
		OPENQASM 3;

		input uint[4] a_in;
		input uint[4] b_in;

		gate majority a, b, c {
			cnot c, b;
			cnot c, a;
			ccnot a, b, c;
		}

		gate unmaj a, b, c {
			ccnot a, b, c;
			cnot c, a;
			cnot a, b;
		}

		qubit cin;
		qubit[4] a;
		qubit[4] b;
		qubit cout;

		// set input states
		for int[8] i in [0: 3] {
		  if(bool(a_in[i])) x a[i];
		  if(bool(b_in[i])) x b[i];
		}

		// add a to b, storing result in b
		majority cin, b[3], a[3];
		for int[8] i in [3: -1: 1] { majority a[i], b[i - 1], a[i - 1]; }
		cnot a[0], cout;
		for int[8] i in [1: 3] { unmaj a[i], b[i - 1], a[i - 1]; }
		unmaj cin, b[3], a[3];

		// todo: subtle bug when trying to get a result type for both at once
		#pragma braket result probability cout, b
		#pragma braket result probability cout
		#pragma braket result probability b
		"""

	grcs_16_qasm = """
		OPENQASM 3.0;

		qubit[16] q;

		h q;
		cz q[0], q[1];
		cz q[6], q[7];
		cz q[8], q[9];
		cz q[14], q[15];

		t q[2:5];
		t q[10:13];

		cz q[4], q[8];
		cz q[6], q[10];

		ry(π/2) q[{0, 1, 14}];
		rx(π/2) q[{7, 9, 15}];

		cz q[1], q[2];
		cz q[9], q[10];

		t q[0];

		ry(π/2) q[4];
		rx(π/2) q[6];

		t q[7];

		rx(π/2) q[8];

		t q[14:15];

		cz q[0], q[4];
		cz q[9], q[13];
		cz q[2], q[6];
		cz q[11], q[15];

		ry(π/2) q[1];

		t q[8];

		ry(π/2) q[10];

		cz q[2], q[3];
		cz q[4], q[5];
		cz q[10], q[11];
		cz q[12], q[13];

		ry(π/2) q[0];

		t q[1];

		rx(π/2) q[6];
		ry(π/2) q[{9, 15}];

		cz q[5], q[9];
		cz q[7], q[11];

		t q[0];

		rx(π/2) q[2];
		ry(π/2) q[3:4];

		t q[6];

		ry(π/2) q[{10, 12}];
		rx(π/2) q[13];

		t q[15];

		cz q[5], q[6];
		cz q[13], q[14];

		t q[2:4];

		ry(π/2) q[{7, 9}];

		t q[10];

		ry(π/2) q[11];

		t q[12];

		cz q[8], q[12];
		cz q[1], q[5];
		cz q[10], q[14];
		cz q[3], q[7];

		rx(π/2) q[6];

		t q[{9, 11}];

		rx(π/2) q[13];

		cz q[0], q[1];
		cz q[6], q[7];
		cz q[8], q[9];
		cz q[14], q[15];

		rx(π/2) q[3];
		ry(π/2) q[{5, 10, 12}];

		t q[13];

		h q;

		#pragma braket result state_vector
		#pragma braket result expectation x(q[0]) 
		#pragma braket result variance x(q[1]) 
	"""

	vqe_qasm = """
		OPENQASM 3.0;
		bit[2] b;
		qubit[2] q;
		ry(-4.97894242803364) q[0];
		ry(-4.435983179322655) q[1];
		cz q[0], q[1];
		ry(6.249142469550989) q[0];
		ry(2.509637558409141) q[1];
		cz q[0], q[1];
		ry(-5.476946260844031) q[0];
		ry(5.937223228770655) q[1];
		cz q[0], q[1];
		ry(1.6839702128823246) q[0];
		ry(-3.211915934619051) q[1];
		b[0] = measure q[0];
		b[1] = measure q[1];
		"""
	all_gates_qasm = """
		OPENQASM 3.0;
		input float theta;
		bit[3] b;
		qubit[3] q;
		rx(0.1) q[0];
		rx(1) q[0];
		prx(0.1, 0.2) q[0];
		x q[0];
		ry(0.1) q[0];
		ry(1) q[0];
		y q[0];
		rz(0.1) q[0];
		rz(1) q[0];
		z q[0];
		h q[0];
		i q[0];
		t q[0];
		ti q[0];
		s q[0];
		si q[0];
		v q[0];
		vi q[0];
		phaseshift(0.1) q[0];
		gpi(0.1) q[0];
		gpi2(0.1) q[0];

		cz q[0], q[1];
		cnot q[0], q[1];
		cy q[0], q[1];
		cv q[0], q[1];

		ecr q[0], q[1];
		swap q[0], q[1];
		iswap q[0], q[1];

		xx(theta) q[0], q[1];
		yy(theta) q[0], q[1];
		xy(theta) q[0], q[1];
		zz(theta) q[0], q[1];
		pswap(theta) q[0], q[1];
		ms(0.1, 0.2, 0.3) q[0], q[1];

		cphaseshift(6.249142469550989) q[0], q[1];
		cphaseshift00(6.249142469550989) q[0], q[1];
		cphaseshift01(6.249142469550989) q[0], q[1];
		cphaseshift10(6.249142469550989) q[0], q[1];

		ccnot q[0], q[1], q[2];
		cswap q[0], q[1], q[2];
		b[0] = measure q[0];
		b[1] = measure q[1];
		b[2] = measure q[2];
		"""
	sv_exact_results_qasm = """
		OPENQASM 3.0;
		qubit[2] q;
		h q;
		#pragma braket result amplitude '00', '01', '10', '11'
		#pragma braket result state_vector
		#pragma braket result density_matrix
		#pragma braket result probability
		#pragma braket result expectation x(q[0])
		#pragma braket result expectation x(q[0]) @ x(q[1])
		#pragma braket result expectation z(q[0]) @ z(q[1])
		#pragma braket result expectation y(q[0]) @ y(q[1])
		#pragma braket result expectation h(q[0]) @ h(q[1])
		#pragma braket result expectation i(q[0]) @ i(q[1])
		#pragma braket result variance x(q[0]) @ y(q[1])
		"""
	dm_exact_results_qasm = """
		OPENQASM 3.0;
		qubit[2] q;
		h q;
		#pragma braket result density_matrix
		#pragma braket result probability
		#pragma braket result expectation x(q[0])
		#pragma braket result variance x(q[0]) @ y(q[1])
		#pragma braket result variance y(q[0])
		"""
	shots_results_qasm = """
		OPENQASM 3.0;
		qubit[10] q;
		h q;
		#pragma braket result probability
		#pragma braket result expectation x(q[0])
		#pragma braket result variance x(q[0]) @ y(q[1])
		#pragma braket result sample x(q[0]) @ y(q[1])
		#pragma braket result expectation z(q[2]) @ z(q[3])
		#pragma braket result expectation x(q[4]) @ x(q[5])
		#pragma braket result expectation y(q[6]) @ y(q[7])
		#pragma braket result expectation h(q[8]) @ h(q[9])
		"""
	@compile_workload begin
		using BraketSimulator, Quasar
		Quasar.builtin_gates[] = BraketSimulator.builtin_gates
		Quasar.parse_pragma[]  = BraketSimulator.parse_pragma
		Quasar.visit_pragma[]  = BraketSimulator.visit_pragma
		simulator              = StateVectorSimulator(5, 0)
		oq3_program            = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), custom_qasm, nothing)
		simulate(simulator, oq3_program, 100)
		simulate("braket_sv_v2", custom_qasm, "{}", 100)

		simulator = DensityMatrixSimulator(2, 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), noise_qasm, nothing)
		simulate(simulator, oq3_program, 100)
		simulate(simulator, [oq3_program, oq3_program], 100)
		simulate("braket_dm_v2", noise_qasm, "{}", 100)
		simulate("braket_dm_v2", [noise_qasm, noise_qasm], "[{}, {}]", 100)

		simulator = StateVectorSimulator(3, 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), unitary_qasm, nothing)
		simulate(simulator, oq3_program, 100)
		simulate(simulator, [oq3_program, oq3_program], 100)
		simulate("braket_sv_v2", unitary_qasm, "{}", 100)
		simulate("braket_sv_v2", [unitary_qasm, unitary_qasm], "[{}, {}]", 100)

		simulator = StateVectorSimulator(6, 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), sv_adder_qasm, Dict("a_in"=>3, "b_in"=>7))
		simulate(simulator, oq3_program, 0)
		simulate("braket_sv_v2", sv_adder_qasm, "{\"a_in\":3,\"b_in\":7}", 100)

		simulator = StateVectorSimulator(16, 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), grcs_16_qasm, nothing)
		simulate(simulator, oq3_program, 0)
		simulate("braket_sv_v2", grcs_16_qasm, "{}", 0)

		simulator = StateVectorSimulator(2, 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), vqe_qasm, nothing)
		simulate(simulator, oq3_program, 100)
		simulate("braket_sv_v2", vqe_qasm, "{}", 100)

		sv_simulator = StateVectorSimulator(3, 0)
		dm_simulator = DensityMatrixSimulator(3, 0)
		oq3_program  = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), all_gates_qasm, Dict("theta"=>0.665))
		simulate(sv_simulator, oq3_program, 100)
		simulate(dm_simulator, oq3_program, 100)
		simulate("braket_sv_v2", all_gates_qasm, "{\"theta\":0.665}", 100)
		simulate("braket_dm_v2", all_gates_qasm, "{\"theta\":0.665}", 100)

		sv_simulator = StateVectorSimulator(2, 0)
		dm_simulator = DensityMatrixSimulator(2, 0)
		sv_oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), sv_exact_results_qasm, nothing)
		dm_oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), dm_exact_results_qasm, nothing)
		simulate(sv_simulator, sv_oq3_program, 0)
		simulate(dm_simulator, dm_oq3_program, 0)
		simulate("braket_sv_v2", sv_exact_results_qasm, "{}", 0)
		simulate("braket_dm_v2", dm_exact_results_qasm, "{}", 0)
		oq3_program = OpenQasmProgram(braketSchemaHeader("braket.ir.openqasm.program", "1"), shots_results_qasm, nothing)
		simulate(sv_simulator, oq3_program, 10)
		simulate(dm_simulator, oq3_program, 10)
		simulate("braket_sv_v2", shots_results_qasm, "{}", 10)
		simulate("braket_dm_v2", shots_results_qasm, "{}", 10)
	end
	"""
		calculate_shots_per_state(branched_sim::BranchedSimulatorOperators) -> Vector{Vector{Int}}

	Takes the shots and instruction sequences from a branched simulator, evolves the states of each instruction,
	and calculates the number of shots allocated to each possible state.

	Returns a list of lists, where each inner list represents a shot and each index of the list represents a qubit.
	For example, for a 3-qubit system, a shot might be represented as [0, 1, 0], indicating qubit 0 is in state 0,
	qubit 1 is in state 1, and qubit 2 is in state 0.
	"""
	function calculate_shots_per_state(branched_sim::BranchedSimulatorOperators)
		# Dictionary to store the shots per state
		shots_per_state = Dict{String, Int}()

		# Process each active path
		for path_idx in branched_sim.active_paths
			# Get the number of shots for this path
			path_shots = branched_sim.shots[path_idx]

			# Skip paths with zero shots
			path_shots <= 0 && continue

			# Calculate the current state by evolving all instructions
			current_state = calculate_current_state(branched_sim, path_idx)

			# Get the number of qubits in the system
			n_qubits = Int(log2(length(current_state)))

			# Calculate probabilities for all possible states
			probs = abs2.(current_state)

			# Create a list of (state_idx, prob) pairs for states with non-zero probability
			state_probs = [(i, p) for (i, p) in enumerate(probs) if p > 0]

			# Sort by probability in descending order
			sort!(state_probs, by = x -> x[2], rev = true)

			# Allocate shots using a more accurate method that preserves the total
			remaining_shots = path_shots
			allocated_shots = Dict{Int, Int}()

			# First pass: allocate integer shots based on probability
			for (state_idx, prob) in state_probs
				if remaining_shots <= 0
					break
				end

				# Calculate shots for this state (floor to ensure we don't over-allocate)
				state_shots = floor(Int, path_shots * prob)

				# Update remaining shots
				remaining_shots -= state_shots

				# Store allocated shots
				allocated_shots[state_idx] = state_shots
			end

			# Second pass: distribute any remaining shots to states with highest fractional part
			if remaining_shots > 0
				# Calculate fractional parts
				fractional_parts = [(i, path_shots * p - floor(path_shots * p)) for (i, p) in state_probs]

				# Sort by fractional part in descending order
				sort!(fractional_parts, by = x -> x[2], rev = true)

				# Distribute remaining shots
				for (state_idx, _) in fractional_parts
					if remaining_shots <= 0
						break
					end

					allocated_shots[state_idx] = get(allocated_shots, state_idx, 0) + 1
					remaining_shots -= 1
				end
			end

			# Convert state indices to binary strings and update shots_per_state
			for (state_idx, state_shots) in allocated_shots
				# Skip states with zero shots
				state_shots <= 0 && continue

				# Convert state index to binary representation (state string)
				# Subtract 1 from state_idx because Julia is 1-indexed but our states start at 0
				state_bits = digits(state_idx - 1, base = 2, pad = n_qubits)
				# Reverse to get the correct endianness
				reverse!(state_bits)
				state_str = join(state_bits)

				# Add the shots for this state to the dictionary
				if haskey(shots_per_state, state_str)
					shots_per_state[state_str] += state_shots
				else
					shots_per_state[state_str] = state_shots
				end
			end
		end

		# Convert the dictionary to the required format: list of lists
		# where each inner list represents a shot and each index represents a qubit
		result = Vector{Vector{Int}}()

		# For each state and its shot count
		for (state_str, shot_count) in shots_per_state
			# Convert the state string to a vector of integers
			state_vec = [parse(Int, c) for c in state_str]

			# Add this state vector to the result shot_count times
			for _ in 1:shot_count
				push!(result, copy(state_vec))
			end
		end

		return result
	end

end # module BraketSimulator
end
