module BraketSimulator

using PrecompileTools: @setup_workload, @compile_workload

using Braket,
    Braket.Observables,
    Dates,
    LinearAlgebra,
    StaticArrays,
    StatsBase,
    Combinatorics,
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

export StateVector, StateVectorSimulator, DensityMatrixSimulator, evolve!, simulate, U, MultiQubitPhaseShift, MultiRZ

const StateVector{T} = Vector{T}
const DensityMatrix{T} = Matrix{T}
const AbstractStateVector{T} = AbstractVector{T}
const AbstractDensityMatrix{T} = AbstractMatrix{T}

abstract type AbstractSimulator <: Braket.AbstractBraketSimulator end
Braket.name(s::AbstractSimulator) = device_id(s)
ap_size(shots::Int, qubit_count::Int) = (shots > 0 && qubit_count < 30) ? 2^qubit_count : 0

include("validation.jl")
include("custom_gates.jl")

const BuiltinGates = merge(Braket.StructTypes.subtypes(Braket.Gate), custom_gates) 

const OBS_LIST = (Observables.X(), Observables.Y(), Observables.Z())
const CHUNK_SIZE = 2^10

parse_program(d, program, shots::Int) = throw(MethodError(parse_program, (d, program, shots)))

function _index_to_endian_bits(ix::Int, qubit_count::Int)
    bits = Vector{Int}(undef, qubit_count)
    for q = 0:qubit_count-1
        b = (ix >> q) & 1
        bits[end-q] = b
    end
    return bits
end

function _formatted_measurements(simulator::D, measured_qubits::Vector{Int}=collect(0:qubit_count(simulator)-1)) where {D<:AbstractSimulator}
    sim_samples = samples(simulator)
    n_qubits    = qubit_count(simulator)
    formatted   = [_index_to_endian_bits(sample, n_qubits)[measured_qubits .+ 1] for sample in sim_samples]
    return formatted
end

function _bundle_results(
    results::Vector{Braket.ResultTypeValue},
    circuit_ir::Program,
    simulator::D;
    measured_qubits::Vector{Int} = collect(0:qubit_count(simulator)-1)
) where {D<:AbstractSimulator}
    task_mtd = Braket.TaskMetadata(
        Braket.header_dict[Braket.TaskMetadata],
        string(uuid4()),
        simulator.shots,
        device_id(simulator),
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
    formatted_samples = simulator.shots > 0 ? _formatted_measurements(simulator, measured_qubits) : Vector{Int}[]
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
    simulator::D;
    measured_qubits::Vector{Int} = collect(0:qubit_count(simulator)-1)
) where {D<:AbstractSimulator}
    task_mtd = Braket.TaskMetadata(
        Braket.header_dict[Braket.TaskMetadata],
        string(uuid4()),
        simulator.shots,
        device_id(simulator),
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
    formatted_samples = simulator.shots > 0 ? _formatted_measurements(simulator, measured_qubits) : Vector{Int}[]
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
    simulator::D,
) where {D<:AbstractSimulator}
    result_values = [calculate(result_type, simulator) for result_type in result_types]
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
    # fix missing targets - what is this? why does it happen?
    for (ri, rr) in enumerate(raw_results)
        if occursin("\"targets\":null", rr)
            raw_results[ri] = replace(
                rr,
                "\"targets\":null" => "\"targets\":$(repr(collect(0:qubit_count-1)))",
            )
        end
    end
    translated_results = [JSON3.read(r, Braket.Result) for r in raw_results]
    return translated_results
end

function _compute_exact_results(d::AbstractSimulator, program::Program, qc::Int, inputs::Dict{String, Float64})
    result_types = _translate_result_types(program.results, qc)
    _validate_result_types_qubits_exist(result_types, qc)
    return _generate_results(program.results, result_types, d)
end

"""
    simulate(simulator::AbstractSimulator, circuit_ir; shots::Int = 0, kwargs...) -> GateModelTaskResult

Simulate the evolution of a statevector or density matrix using the passed in `simulator`.
The instructions to apply (gates and noise channels) and measurements to make are
encoded in `circuit_ir`. Supported IR formats are `OpenQASMProgram` (OpenQASM3)
and `Program` (JAQCD). Returns a `GateModelTaskResult` containing the individual shot
measurements (if `shots > 0`), final calculated results, circuit IR, and metadata
about the task.
"""
function simulate(
    simulator::AbstractSimulator,
    circuit_ir::OpenQasmProgram;
    shots::Int = 0,
    kwargs...,
)
    program = parse_program(simulator, circuit_ir, shots)
    qc      = qubit_count(program)
    _validate_ir_results_compatibility(simulator, program.results, Val(:JAQCD))
    _validate_ir_instructions_compatibility(simulator, program, Val(:JAQCD))
    _validate_shots_and_ir_results(shots, program.results, qc)
    operations = program.instructions
    if shots > 0 && !isempty(program.basis_rotation_instructions)
        operations = vcat(operations, program.basis_rotation_instructions)
    end
    inputs        = isnothing(circuit_ir.inputs) ? Dict{String, Float64}() : Dict{String, Float64}(k=>v for (k,v) in circuit_ir.inputs)
    symbol_inputs = Dict{Symbol,Number}(Symbol(k) => v for (k, v) in inputs)
    operations    = [Braket.bind_value!(op, symbol_inputs) for op in operations]
    _validate_operation_qubits(operations)
    reinit!(simulator, qc, shots)
    stats = @timed begin
        simulator = evolve!(simulator, operations)
    end
    @debug "Time for evolution: $(stats.time)"
    results = shots == 0 && !isempty(program.results) ? _compute_exact_results(simulator, program, qc, inputs) : [Braket.ResultTypeValue(rt, 0.0) for rt in program.results]
    mqs     = get(kwargs, :measured_qubits, collect(0:qc-1))
    isempty(mqs) && (mqs = collect(0:qc-1))
    stats   = @timed _bundle_results(results, circuit_ir, simulator; measured_qubits=mqs)
    @debug "Time for results bundling: $(stats.time)"
    res = stats.value
    return res
end

function simulate(
    simulator::AbstractSimulator,
    circuit_ir::Program,
    qubit_count::Int;
    shots::Int = 0,
    kwargs...,
)
    _validate_ir_results_compatibility(simulator, circuit_ir.results, Val(:JAQCD))
    _validate_ir_instructions_compatibility(simulator, circuit_ir, Val(:JAQCD))
    _validate_shots_and_ir_results(shots, circuit_ir.results, qubit_count)
    operations = circuit_ir.instructions
    if shots > 0 && !isempty(circuit_ir.basis_rotation_instructions)
        operations = vcat(operations, circuit_ir.basis_rotation_instructions)
    end
    inputs = get(kwargs, :inputs, Dict{String,Float64}())
    symbol_inputs = Dict{Symbol,Number}(Symbol(k) => v for (k, v) in inputs)
    operations = [bind_value!(Instruction(op), symbol_inputs) for op in operations]
    _validate_operation_qubits(operations)
    reinit!(simulator, qubit_count, shots)
    stats = @timed begin
        simulator = evolve!(simulator, operations)
    end
    @debug "Time for evolution: $(stats.time)"
    results = shots == 0 && !isempty(circuit_ir.results) ? _compute_exact_results(simulator, circuit_ir, qubit_count, inputs) : Braket.ResultTypeValue[]
    mqs     = get(kwargs, :measured_qubits, collect(0:qubit_count-1))
    isempty(mqs) && (mqs = collect(0:qubit_count-1))
    stats = @timed _bundle_results(results, circuit_ir, simulator; measured_qubits=mqs)
    @debug "Time for results bundling: $(stats.time)"
    res = stats.value
    return res
end

include("gate_kernels.jl")
include("noise_kernels.jl")
include("observables.jl")
include("result_types.jl")
include("properties.jl")
include("inverted_gates.jl")
include("pow_gates.jl")
include("sv_simulator.jl")
include("dm_simulator.jl")

function __init__()
    Braket._simulator_devices[]["braket_dm_v2"] =
        DensityMatrixSimulator{ComplexF64,DensityMatrix{ComplexF64}}
    Braket._simulator_devices[]["braket_sv_v2"] =
        StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
    # Braket._simulator_devices[]["default"] =
    #    StateVectorSimulator{ComplexF64,StateVector{ComplexF64}}
end

end # module BraketSimulator
