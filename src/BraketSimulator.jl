module BraketSimulator

using PrecompileTools: @setup_workload, @compile_workload

using Braket,
    Braket.Observables,
    Dates,
    LinearAlgebra,
    StaticArrays,
    StatsBase,
    Combinatorics,
    OpenQASM3,
    UUIDs,
    JSON3,
    Random

import Braket:
    Instruction,
    AbstractBraketSimulator,
    Program,
    OpenQasmProgram,
    ir_typ,
    apply_gate!,
    apply_noise!,
    qubit_count,
    I,
    simulate,
    device_id,
    bind_value!

export StateVector, StateVectorSimulator, DensityMatrixSimulator, evolve!, classical_shadow, simulate

const StateVector{T} = Vector{T}
const DensityMatrix{T} = Matrix{T}
const AbstractStateVector{T} = AbstractVector{T}
const AbstractDensityMatrix{T} = AbstractMatrix{T}

abstract type AbstractSimulator <: Braket.AbstractBraketSimulator end
Braket.name(s::AbstractSimulator) = device_id(s)
ap_size(shots::Int, qubit_count::Int) = (shots > 0 && qubit_count < 30) ? 2^qubit_count : 0

include("validation.jl")
include("builtins.jl")
include("custom_gates.jl")
include("openqasm.jl")

const BuiltinGates = merge(Braket.StructTypes.subtypes(Braket.Gate), custom_gates) 

const OBS_LIST = (Observables.X(), Observables.Y(), Observables.Z())
const CHUNK_SIZE = 2^10

function parse_program(d::D, program::OpenQasmProgram, shots::Int) where {D<:AbstractSimulator}
    if endswith(program.source, "qasm") && isfile(program.source)
        parsed_prog      = OpenQASM3.parse(read(program.source, String))
    else
        parsed_prog      = OpenQASM3.parse(program.source)
    end
    interpreted_circ = interpret(parsed_prog, program.inputs)
    if shots > 0
        Braket.validate_circuit_and_shots(interpreted_circ, shots)
        Braket.basis_rotation_instructions!(interpreted_circ)
    end
    return convert(Braket.Program, interpreted_circ)
end

function index_to_endian_bits(ix::Int, qc::Int)
    bits = Vector{Int}(undef, qc)
    for q = 0:qc-1
        b = (ix >> q) & 1
        bits[end-q] = b
    end
    return bits
end

function _formatted_measurements(d::D) where {D<:AbstractSimulator}
    sim_samples = samples(d)
    qc          = qubit_count(d)
    formatted   = [index_to_endian_bits(sample, qc) for sample in sim_samples]
    return formatted
end

_get_measured_qubits(qc::Int) = collect(0:qc-1)
function _bundle_results(
    results::Vector{Braket.ResultTypeValue},
    circuit_ir::Program,
    d::D,
) where {D<:AbstractSimulator}
    task_mtd = Braket.TaskMetadata(
        Braket.header_dict[Braket.TaskMetadata],
        string(uuid4()),
        d.shots,
        device_id(d),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    addl_mtd = Braket.AdditionalMetadata(
        circuit_ir,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    formatted_samples = d.shots > 0 ? _formatted_measurements(d) : Vector{Int}[]
    measured_qubits   = collect(0:qubit_count(d)-1)
    return Braket.GateModelTaskResult(
        Braket.header_dict[Braket.GateModelTaskResult],
        formatted_samples,
        nothing,
        results,
        measured_qubits,
        task_mtd,
        addl_mtd,
    )
end

function _bundle_results(
    results::Vector{Braket.ResultTypeValue},
    circuit_ir::OpenQasmProgram,
    d::D,
) where {D<:AbstractSimulator}
    task_mtd = Braket.TaskMetadata(
        Braket.header_dict[Braket.TaskMetadata],
        string(uuid4()),
        d.shots,
        device_id(d),
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    addl_mtd = Braket.AdditionalMetadata(
        circuit_ir,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
        nothing,
    )
    formatted_samples = d.shots > 0 ? _formatted_measurements(d) : Vector{Int}[]
    measured_qubits   = collect(0:qubit_count(d)-1)
    return Braket.GateModelTaskResult(
        Braket.header_dict[Braket.GateModelTaskResult],
        formatted_samples,
        nothing,
        results,
        measured_qubits,
        task_mtd,
        addl_mtd,
    )
end

function _generate_results(
    results::Vector{<:Braket.AbstractProgramResult},
    result_types::Vector,
    d::D,
) where {D<:AbstractSimulator}
    result_values = [calculate(result_type, d) for result_type in result_types]
    result_values =
        [val isa Matrix ? Braket.complex_matrix_to_ir(val) : val for val in result_values]
    return [
        Braket.ResultTypeValue(result, result_value) for
        (result, result_value) in zip(results, result_values)
    ]
end

function _translate_result_types(
    results::Vector{Braket.AbstractProgramResult},
    qubit_count::Int,
)
    # results are in IR format
    raw_results = [JSON3.write(r) for r in results]
    # fix missing targets
    for (ri, rr) in enumerate(raw_results)
        if occursin("\"targets\":null", rr)
            raw_results[ri] = replace(
                rr,
                "\"targets\":null" => "\"targets\":$(repr(collect(0:qubit_count-1)))",
            )
        end
    end
    is_ag = results[1] isa Braket.IR.AdjointGradient
    translated_results =
        is_ag ? [JSON3.read(r, Braket.AdjointGradient) for r in raw_results] :
        [JSON3.read(r, Braket.Result) for r in raw_results]
    return translated_results
end

function classical_shadow(
    d::LocalSimulator,
    obs_qubits::Vector{Int},
    circuit::Program,
    shots::Int,
    seed::Int,
)
    n_qubits = length(obs_qubits)
    n_snapshots = shots
    recipes = rand(Xoshiro(seed), 0:2, (n_snapshots, n_qubits))
    outcomes = compute_shadows(d._delegate, circuit, recipes, obs_qubits)
    return outcomes, recipes
end
function classical_shadow(
    d::LocalSimulator,
    obs_qubits::Vector{Vector{Int}},
    circuits::Vector{Program},
    shots::Int,
    seed::Vector{Int},
)
    all_outcomes =
        [zeros(Float64, (shots, length(obs_qubits[ix]))) for ix = 1:length(circuits)]
    all_recipes = [zeros(Int, (shots, length(obs_qubits[ix]))) for ix = 1:length(circuits)]
    Threads.@threads for ix = 1:length(circuits)
        start = time()
        outcomes, recipes =
            classical_shadow(d, obs_qubits[ix], circuits[ix], shots, seed[ix])
        stop = time()
        @views begin
            all_outcomes[ix] .= outcomes
            all_recipes[ix] .= recipes
        end
    end
    return all_outcomes, all_recipes
end

function compute_shadows(
    d::AbstractSimulator,
    circuit_ir::Program,
    recipes,
    obs_qubits;
    kwargs...,
)
    qc = qubit_count(circuit_ir)
    _validate_ir_instructions_compatibility(d, circuit_ir, Val(:JAQCD))
    operations = circuit_ir.instructions
    operations = [Instruction(op) for op in operations]
    _validate_operation_qubits(operations)

    n_snapshots = size(recipes, 1)
    n_qubits = size(recipes, 2)
    reinit!(d, qc, n_snapshots)
    d = evolve!(d, operations)
    snapshot_d = similar(d, shots = 1)
    qubit_inds =
        [findfirst(q -> q == oq, 0:maximum(qubits(circuit_ir))) for oq in obs_qubits]
    any(isnothing, qubit_inds) && throw(
        ErrorException(
            "one of obs_qubits ($obs_qubits) not present in circuit with qubits $(sort(collect(qubits(circuit_ir)))).",
        ),
    )
    measurements = Matrix{Int}(undef, n_snapshots, n_qubits)
    for s_ix = 1:n_snapshots
        copyto!(snapshot_d, d)
        snapshot_rotations = reduce(
            vcat,
            [
                diagonalizing_gates(OBS_LIST[recipes[s_ix, wire_idx]+1], [wire]) for
                (wire_idx, wire) in enumerate(obs_qubits)
            ],
        )
        snapshot_d = evolve!(snapshot_d, snapshot_rotations)
        @views begin
            measurements[s_ix, :] .= _formatted_measurements(snapshot_d)[1][qubit_inds]
        end
    end
    return measurements
end

function _compute_exact_results(d::AbstractSimulator, program::Program, qc::Int, inputs::Dict{String, Float64})
    result_types = _translate_result_types(program.results, qc)
    _validate_result_types_qubits_exist(result_types, qc)
    if program.results[1] isa Braket.IR.AdjointGradient
        rt = result_types[1]
        ev, grad =
            calculate(rt, Vector{Instruction}(program.instructions), inputs, d)
        result_val = Dict{String,Any}("expectation" => ev, "gradient" => grad)
        return [Braket.ResultTypeValue(program.results[1], result_val)]
    else
        return _generate_results(program.results, result_types, d)
    end
end

function simulate(
    d::AbstractSimulator,
    circuit_ir::OpenQasmProgram;
    shots::Int = 0,
    kwargs...,
)
    program = parse_program(d, circuit_ir, shots)
    qc      = qubit_count(program)
    _validate_ir_results_compatibility(d, program.results, Val(:JAQCD))
    _validate_ir_instructions_compatibility(d, program, Val(:JAQCD))
    _validate_shots_and_ir_results(shots, program.results, qc)
    operations = program.instructions
    if shots > 0 && !isempty(program.basis_rotation_instructions)
        operations = vcat(operations, program.basis_rotation_instructions)
    end
    inputs        = isnothing(circuit_ir.inputs) ? Dict{String, Float64}() : Dict{String, Float64}(k=>v for (k,v) in circuit_ir.inputs)
    symbol_inputs = Dict{Symbol,Number}(Symbol(k) => v for (k, v) in inputs)
    operations    = [Braket.bind_value!(op, symbol_inputs) for op in operations]
    _validate_operation_qubits(operations)
    reinit!(d, qc, shots)
    stats = @timed begin
        d = evolve!(d, operations)
    end
    @debug "Time for evolution: $(stats.time)"
    results = shots == 0 && !isempty(program.results) ? _compute_exact_results(d, program, qc, inputs) : [Braket.ResultTypeValue(rt, 0.0) for rt in program.results]
    stats   = @timed _bundle_results(results, circuit_ir, d)
    @debug "Time for results bundling: $(stats.time)"
    res = stats.value
    return res
end

function simulate(
    d::AbstractSimulator,
    circuit_ir::Program,
    qc::Int;
    shots::Int = 0,
    kwargs...,
)
    _validate_ir_results_compatibility(d, circuit_ir.results, Val(:JAQCD))
    _validate_ir_instructions_compatibility(d, circuit_ir, Val(:JAQCD))
    _validate_shots_and_ir_results(shots, circuit_ir.results, qc)
    operations = circuit_ir.instructions
    if shots > 0 && !isempty(circuit_ir.basis_rotation_instructions)
        operations = vcat(operations, circuit_ir.basis_rotation_instructions)
    end
    inputs = get(kwargs, :inputs, Dict{String,Float64}())
    symbol_inputs = Dict{Symbol,Number}(Symbol(k) => v for (k, v) in inputs)
    operations = [bind_value!(Instruction(op), symbol_inputs) for op in operations]
    _validate_operation_qubits(operations)
    reinit!(d, qc, shots)
    stats = @timed begin
        d = evolve!(d, operations)
    end
    @debug "Time for evolution: $(stats.time)"
    results = shots == 0 && !isempty(circuit_ir.results) ? _compute_exact_results(d, circuit_ir, qc, inputs) : Braket.ResultTypeValue[]
    stats = @timed _bundle_results(results, circuit_ir, d)
    @debug "Time for results bundling: $(stats.time)"
    res = stats.value
    return res
end

include("gate_kernels.jl")
include("noise_kernels.jl")
include("observables.jl")
include("result_types.jl")
include("properties.jl")
include("derivative_gates.jl")
include("inverted_gates.jl")
include("pow_gates.jl")
include("sv_simulator.jl")
include("dm_simulator.jl")

function __init__()
    Braket._simulator_devices[]["braket_jl_dm"] =
        DensityMatrixSimulator{ComplexF64,DensityMatrix{ComplexF64}}
    Braket._simulator_devices[]["braket_jl_sv"] =
        StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
    Braket._simulator_devices[]["default"] =
        StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
end


@setup_workload begin
    @compile_workload begin
        function ghz_circuit(qubit_count::Int)
           ghz_circ = Braket.Circuit()
           Braket.H(ghz_circ, 0)
           for target_qubit = 1:qubit_count-1
               Braket.CNot(ghz_circ, 0, target_qubit)
           end
           return ghz_circ
        end
        n_qubits = 10
        svs = StateVectorSimulator(n_qubits, 0)
        oq3_src = "OPENQASM 3.0;\nbit[10] b;\nqubit[10] q;\nh q[0];\ncnot q[0], q[1];\ncnot q[0], q[2];\ncnot q[0], q[3];\ncnot q[0], q[4];\ncnot q[0], q[5];\ncnot q[0], q[6];\ncnot q[0], q[7];\ncnot q[0], q[8];\ncnot q[0], q[9];\nb[0] = measure q[0];\nb[1] = measure q[1];\nb[2] = measure q[2];\nb[3] = measure q[3];\nb[4] = measure q[4];\nb[5] = measure q[5];\nb[6] = measure q[6];\nb[7] = measure q[7];\nb[8] = measure q[8];\nb[9] = measure q[9];"
        simulate(svs, Braket.OpenQasmProgram(Braket.header_dict[Braket.OpenQasmProgram], oq3_src, nothing); shots=10)
    end
end
#include("precompile.jl")

end # module BraketCircuitSimulator
