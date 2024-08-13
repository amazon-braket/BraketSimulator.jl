module BraketSimulatorPythonExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, PythonCall, BraketSimulator.Dates
end

import BraketSimulator:
    qubit_count,
    Instruction,
    Observables,
    AbstractProgramResult,
    ResultTypeValue,
    GateModelTaskResult,
    Program,
    Gate,
    AngledGate,
    AbstractIR,
    AbstractProgram,
    IRObservable,
    AbstractSimulator,
    simulate,
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
    PythonCall.GIL.@lock begin
        braket[]    = pyimport("braket")
        numpy[]     = pyimport("numpy")
        pyimport("braket.ir.openqasm.program_v1")
        pyimport("braket.task_result.task_metadata_v1")
        pyimport("braket.task_result.additional_metadata")
    end
end

@recompile_invalidations begin
    function simulate(simulator::AbstractSimulator,
                      task_specs::PyList{Any},
                      args...;
                      max_parallel::Int=-1,
                      inputs::Union{PyDict{Any, Any}, PyList{Any}, Nothing} = nothing,
                      shots::Int=0,
                      kwargs...
                     )
        jl_inputs = Dict{String, Float64}()
        jl_specs = Vector{Union{BraketSimulator.OpenQasmProgram, BraketSimulator.Program}}(undef, length(task_specs))
        PythonCall.GIL.@lock begin
            jl_specs  = convert(Vector{Union{BraketSimulator.OpenQasmProgram, BraketSimulator.Program}}, task_specs)
            if !isnothing(inputs)
                py_inputs = inputs isa PyDict ? [inputs] : inputs
                n_inputs  = length(py_inputs)
                jl_inputs = Vector{Dict{String, Float64}}(undef, n_inputs)
                for input_ix in 1:n_inputs
                    jl_inputs[input_ix] = Dict{String, Float64}(pyconvert(String, k)=>pyconvert(Float64, v) for (k,v) in py_inputs[input_ix])
                end
            end
        end
        jl_kwargs = (max_parallel=max_parallel, inputs=jl_inputs)
        jl_results = simulate(simulator, jl_specs, shots; jl_kwargs...)
        return PythonCall.GIL.@lock map(Py, jl_results)
    end
end

end
