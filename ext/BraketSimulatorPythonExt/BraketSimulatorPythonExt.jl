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
                      task_args::Union{PyList{Any}, Nothing},
                      task_kwargs::Union{PyList{Any}, Nothing};
                      max_parallel::Int=-1,
                     )
        args_provided = !isnothing(task_args) && !isempty(task_args)
        kwargs_provided = !isnothing(task_kwargs) && !isempty(task_kwargs)
        kwargs_provided && length(task_kwargs) > 1 && length(task_kwargs) != length(task_specs) && throw(ArgumentError("number of kwargs ($(length(task_kwargs))) is not equal to number of task specifications ($(length(task_specs))) or equal to 1.")) 
        args_provided && length(task_args) > 1 && length(task_args) != length(task_specs) && throw(ArgumentError("number of args ($(length(task_args))) is not equal to number of task specifications ($(length(task_specs))) or equal to 1.")) 
        PythonCall.GC.disable()
        jl_specs = convert(Vector{Union{Braket.OpenQasmProgram, Braket.Program}}, task_specs)
        jl_args  = args_provided ? map(arg->pyconvert(Vector{Any}, arg), task_args) : [[] for t_ix in 1:length(task_specs)]
        
        inputs = Dict{String, Float64}()
        if kwargs_provided
            inputs = Vector{Dict{String, Float64}}(undef, length(task_specs))
            # we have to extract `shots` from the kwargs lists if present
            for (task_ix, task_kwarg) in enumerate(task_kwargs)
                py_inputs = get(task_kwarg, "inputs", [])
                inputs[task_ix] = Dict{String, Float64}(pyconvert(String, k)=>pyconvert(Float64, v) for (k,v) in py_inputs)
                if pyhasitem(task_kwarg, "shots")
                    jl_shots = pyconvert(Int, task_kwarg["shots"])
                    isempty(jl_args[task_ix]) && push!(jl_args[task_ix], jl_shots)
                    jl_args[task_ix][end] = jl_shots
                end
            end
        end
        jl_results = simulate(simulator, jl_specs, jl_args; max_parallel=max_parallel, inputs=inputs)
        PythonCall.GC.enable()
        return map(Py, jl_results)
    end
end

end
