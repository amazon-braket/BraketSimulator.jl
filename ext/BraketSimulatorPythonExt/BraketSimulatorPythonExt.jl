module BraketSimulatorPythonExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, PythonCall, JSON3
end

import BraketSimulator: simulate

function simulate(simulator, task_spec::PyList{Any}, shots; inputs=Dict{String,Float64}(), kwargs...)
    jl_inputs  = inputs isa Dict ? inputs : [pyconvert(Dict{String, Any}, input) for input in inputs]
    jl_results = simulate(simulator, BraketSimulator.OpenQasmProgram[t for t in task_spec], shots; inputs=jl_inputs, kwargs...)
    jsons      = (JSON3.write(r) for r in jl_results)
    py_results = pylist(jsons)
    return py_results
end

end
