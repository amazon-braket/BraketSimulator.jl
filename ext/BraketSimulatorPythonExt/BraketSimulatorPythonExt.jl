module BraketSimulatorPythonExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, PythonCall, JSON3
end

using BraketSimulator: simulate, AbstractProgramResult

function BraketSimulator.simulate(simulator_id::String, task_spec::String, py_inputs::String, shots::Int; kwargs...)
    inputs    = JSON3.read(py_inputs, Dict{String, Any})
    jl_spec   = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, inputs)
    simulator = if simulator_id == "braket_sv_v2"
        StateVectorSimulator(0, shots)
    elseif simulator_id == "braket_dm_v2"
        DensityMatrixSimulator(0, shots)
    end
    jl_results = @time "simulate" simulate(simulator, jl_spec, shots; kwargs...)
    py_results = @time "JSON write" JSON3.write(jl_results)
    # this is expensive due to allocations
    simulator  = nothing
    inputs     = nothing
    jl_spec    = nothing
    jl_results = nothing
    return py_results
end
function BraketSimulator.simulate(simulator_id::String, task_specs::PyList, py_inputs::String, shots::Int; kwargs...)
    inputs    = JSON3.read(py_inputs, Vector{Dict{String, Any}})
    jl_specs  = map(zip(task_specs, inputs)) do (task_spec, input)
        BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), task_spec, input)
    end
    simulator = if simulator_id == "braket_sv_v2"
        StateVectorSimulator(0, shots)
    elseif simulator_id == "braket_dm_v2"
        DensityMatrixSimulator(0, shots)
    end
    jl_results = simulate(simulator, jl_specs, shots; kwargs...)
    jsons      = [JSON3.write(r) for r in jl_results]
    simulator  = nothing
    jl_results = nothing
    inputs     = nothing
    jl_specs   = nothing
    return jsons
end

end
