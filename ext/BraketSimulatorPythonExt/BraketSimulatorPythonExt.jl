module BraketSimulatorPythonExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, PythonCall, JSON3
end

using BraketSimulator: simulate

function BraketSimulator.simulate(simulator, task_spec::String, inputs::Dict{String, Any}, shots::Int; kwargs...)
    jl_specs   = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, inputs)
    jl_results = simulate(simulator, jl_specs, shots; kwargs...)
    json       = JSON3.write(jl_results)
    return json
end
function BraketSimulator.simulate(simulator, task_specs::Vector{String}, inputs::Vector{Dict{String, Any}}, shots::Int; kwargs...)
    jl_specs   = map(zip(task_specs, inputs)) do (task_spec, input)
        BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, input)
    end
    jl_results = simulate(simulator, jl_specs, shots; kwargs...)
    jsons      = [JSON3.write(r) for r in jl_results]
    return jsons
end

end
