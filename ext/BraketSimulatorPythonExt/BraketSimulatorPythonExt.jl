module BraketSimulatorPythonExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, BraketSimulator.Braket, PythonCall, BraketSimulator.Dates
end

import BraketSimulator.Braket:
    LocalSimulator,
    qubit_count,
    Instruction,
    Observables,
    AbstractProgramResult,
    ResultTypeValue,
    format_result,
    GateModelQuantumTaskResult,
    GateModelTaskResult,
    Program,
    Gate,
    AngledGate,
    AbstractIR,
    AbstractProgram,
    IRObservable
import BraketSimulator:
    AbstractSimulator,
    simulate,
    parse_program,
    DoubleExcitation,
    SingleExcitation,
    Control,
    MultiQubitPhaseShift,
    MultiRZ

const numpy     = Ref{Py}()
const braket    = Ref{Py}()
const sympy     = Ref{Py}()

include("translation.jl")

function __init__()
    # must set these when this code is actually loaded
    braket[]    = pyimport("braket")
    numpy[]     = pyimport("numpy")
    pyimport("braket.ir.openqasm.program_v1")
    pyimport("braket.task_result.task_metadata_v1")
    pyimport("braket.task_result.additional_metadata")
end

@recompile_invalidations begin
    function simulate(simulator::AbstractSimulator,
             task_specs::PyList{Any},
             shots::Int;
             max_parallel::Int=-1,
             inputs::Union{Vector{Dict{String, Float64}}, Dict{String, Float64}} = Dict{String, Float64}(),
             kwargs...
            )
        jl_results = simulate(simulator, convert(Vector{Union{Braket.OpenQasmProgram, Braket.Program}}, task_specs), shots; max_parallel=max_parallel, inputs=inputs, kwargs...)
        return map(Py, jl_results)
    end
end

end
