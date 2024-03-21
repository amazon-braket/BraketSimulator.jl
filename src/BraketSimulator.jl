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
    apply_gate!,
    apply_noise!,
    qubit_count,
    I,
    device_id,
    bind_value!

export StateVector, StateVectorSimulator, DensityMatrixSimulator, evolve!, classical_shadow

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
    interpreted_circ = interpret(parsed_prog, extern_lookup=program.inputs)
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

function (d::AbstractSimulator)(
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
    operations    = [bind_value!(Instruction(op), symbol_inputs) for op in operations]
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

function (d::AbstractSimulator)(
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
    Braket._simulator_devices[]["braket_dm"] =
        DensityMatrixSimulator{ComplexF64,DensityMatrix{ComplexF64}}
    Braket._simulator_devices[]["braket_sv"] =
        StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
    Braket._simulator_devices[]["default"] =
        StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
end

@setup_workload begin
    @compile_workload begin
        for (instructions, qubit_count, state_vector, probability_amplitudes) in [
            ([Instruction(H(), [0])], 1, [0.70710678, 0.70710678], [0.5, 0.5]),
            ([Instruction(X(), [0])], 1, [0, 1], [0, 1]),
            ([Instruction(X(), [0])], 2, [0, 0, 1, 0], [0, 0, 1, 0]),
            ([Instruction(Y(), [0])], 1, [0, im], [0, 1]),
            ([Instruction(X(), [0]), Instruction(X(), [1])], 2, [0, 0, 0, 1], [0, 0, 0, 1]),
            ([Instruction(X(), [0]), Instruction(Z(), [0])], 1, [0, -1], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(CNot(), [0, 1])],
                2,
                [0, 0, 0, 1],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(CY(), [0, 1])],
                2,
                [0, 0, 0, im],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(CZ(), [0, 1])],
                2,
                [0, 0, 1, 0],
                [0, 0, 1, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(Swap(), [0, 1])],
                2,
                [0, 1, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(ISwap(), [0, 1])],
                2,
                [0, im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(Swap(), [0, 2])],
                3,
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
            ),
            ([Instruction(X(), [0]), Instruction(S(), [0])], 1, [0, im], [0, 1]),
            ([Instruction(X(), [0]), Instruction(Si(), [0])], 1, [0, -im], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(T(), [0])],
                1,
                [0, 0.70710678 + 0.70710678 * im],
                [0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(Ti(), [0])],
                1,
                [0, 0.70710678 - 0.70710678 * im],
                [0, 1],
            ),
            ([Instruction(V(), [0])], 1, [0.5 + 0.5 * im, 0.5 - 0.5 * im], [0.5, 0.5]),
            ([Instruction(Vi(), [0])], 1, [0.5 - 0.5 * im, 0.5 + 0.5 * im], [0.5, 0.5]),
            ([Instruction(I(), [0])], 1, [1, 0], [1, 0]),
            ([Instruction(Unitary([0 1; 1 0]), [0])], 1, [0, 1], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(PhaseShift(0.15), [0])],
                1,
                [0, 0.98877108 + 0.14943813 * im],
                [0, 1],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CPhaseShift(0.15), [0, 1]),
                ],
                2,
                [0, 0, 0, 0.98877108 + 0.14943813 * im],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(CPhaseShift00(0.15), [0, 1])],
                2,
                [0.98877108 + 0.14943813 * im, 0, 0, 0],
                [1, 0, 0, 0],
            ),
            (
                [Instruction(X(), [1]), Instruction(CPhaseShift01(0.15), [0, 1])],
                2,
                [0, 0.98877108 + 0.14943813 * im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(CPhaseShift10(0.15), [0, 1])],
                2,
                [0, 0, 0.98877108 + 0.14943813 * im, 0],
                [0, 0, 1, 0],
            ),
            (
                [Instruction(Rx(0.15), [0])],
                1,
                [0.99718882, -0.07492971 * im],
                [0.99438554, 0.00561446],
            ),
            (
                [Instruction(X(), [0]), Instruction(Ry(0.15), [0])],
                1,
                [-0.07492971, 0.99718882],
                [0.00561446, 0.99438554],
            ),
            (
                [Instruction(H(), [0]), Instruction(Rz(0.15), [0])],
                1,
                [0.70511898 - 0.0529833 * im, 0.70511898 + 0.0529833 * im],
                [0.5, 0.5],
            ),
            (
                [Instruction(X(), [0]), Instruction(PSwap(0.15), [0, 1])],
                2,
                [0, 0.98877108 + 0.14943813 * im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(XY(0.15), [0, 1])],
                2,
                [0, 0.07492971 * im, 0.99718882, 0],
                [0, 0.00561446, 0.99438554, 0],
            ),
            (
                [Instruction(XX(0.3), [0, 1])],
                2,
                [0.98877108, 0, 0, -0.14943813 * im],
                [0.97766824, 0, 0, 0.02233176],
            ),
            (
                [Instruction(YY(0.3), [0, 1])],
                2,
                [0.98877108, 0, 0, 0.14943813 * im],
                [0.97766824, 0, 0, 0.02233176],
            ),
            (
                [Instruction(ZZ(0.15), [0, 1])],
                2,
                [0.99718882 - 0.07492971 * im, 0, 0, 0],
                [1, 0, 0, 0],
            ),
            (
                [Instruction(BraketSimulator.MultiRZ(0.15), [0, 1, 2])],
                3,
                [0.99718882 - 0.07492971 * im, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CCNot(), [0, 1, 2]),
                ],
                3,
                [0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 0, 1],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CSwap(), [0, 1, 2]),
                ],
                3,
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
            ),
        ]
            simulation = StateVectorSimulator(qubit_count, 0)
            simulation = evolve!(simulation, instructions)
        end
        for (obs, equivalent_gates, qubit_count) in [
            ([(Braket.Observables.X(), [0])], [Instruction(H(), [0])], 1),
            ([(Braket.Observables.Z(), [0])], Instruction[], 1),
            ([(Braket.Observables.I(), [0])], Instruction[], 1),
            (
                [
                    (Braket.Observables.X(), [0]),
                    (Braket.Observables.Z(), [3]),
                    (Braket.Observables.H(), [2]),
                ],
                [Instruction(H(), [0]), Instruction(Ry(-π / 4), [2])],
                5,
            ),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.X(),
                        Braket.Observables.Z(),
                        Braket.Observables.H(),
                        Braket.Observables.I(),
                    ]),
                    (0, 3, 2, 1),
                )],
                [Instruction(H(), [0]), Instruction(Ry(-π / 4), [2])],
                5,
            ),
            (
                [(Braket.Observables.X(), [0, 1])],
                [Instruction(H(), [0]), Instruction(H(), [1])],
                2,
            ),
            ([(Braket.Observables.Z(), [0, 1])], Instruction[], 2),
            ([(Braket.Observables.I(), [0, 1])], Instruction[], 2),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.I(),
                        Braket.Observables.Z(),
                    ]),
                    (2, 0),
                )],
                Instruction[],
                3,
            ),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.X(),
                        Braket.Observables.Z(),
                    ]),
                    (2, 0),
                )],
                [Instruction(H(), [2])],
                3,
            ),
        ]
            sim_observables = StateVectorSimulator(qubit_count, 0)
            sim_observables = BraketSimulator.apply_observables!(sim_observables, obs)
            sim_gates = StateVectorSimulator(qubit_count, 0)
            sim_gates = BraketSimulator.evolve!(sim_gates, equivalent_gates)
        end

        function qft_circuit_operations(qubit_count::Int)
            qft_ops = Instruction[]
            for target_qubit = 0:qubit_count-1
                angle = π / 2
                push!(qft_ops, Instruction(H(), [target_qubit]))
                for control_qubit = target_qubit+1:qubit_count-1
                    push!(
                        qft_ops,
                        Instruction(CPhaseShift(angle), [control_qubit, target_qubit]),
                    )
                    angle /= 2
                end
            end
            return qft_ops
        end

        qubit_count = 16
        simulation = StateVectorSimulator(qubit_count, 0)
        operations = qft_circuit_operations(qubit_count)
        simulation = BraketSimulator.evolve!(simulation, operations)
        simulation = StateVectorSimulator(2, 10000)
        simulation = BraketSimulator.evolve!(
            simulation,
            [Instruction(H(), [0]), Instruction(CNot(), [0, 1])],
        )
        samples = BraketSimulator.samples(simulation)

        for N = 1:16
            pe = Braket.PauliEigenvalues(Val(N))
            for i = 1:2^N
                pe[i]
            end
        end
    end
end


end # module BraketCircuitSimulator
