convert_targets(::Nothing)       = PythonCall.pybuiltins.None
convert_targets(ts::Vector{Int}) = pylist(ts)

function inputs_to_jl(x)
    if pyis(x, pybuiltins.None)
        return nothing
    else
        jl_inputs = map(x.items()) do x_
            k, v = x_
            jl_k = pyconvert(String, k)
            jl_v = pyisinstance(v, pybuiltins.int) ? pyconvert(Int, v) : pyconvert(Float64, v)
            return jl_k=>jl_v
        end
        return Dict(jl_inputs)
    end
end
function jl_convert_oqprogram(t::Type{OpenQasmProgram}, x::Py)
    bsh    = pyconvert(Braket.braketSchemaHeader, x.braketSchemaHeader)
    source = pyconvert(String, x.source)
    inputs = inputs_to_jl(x.inputs) 
    PythonCall.pyconvert_return(t(bsh, source, inputs))
end

function jl_convert_bsh(t::Type{Braket.braketSchemaHeader}, x::Py)
    name    = pyconvert(String, x.name)
    version = pyconvert(String, x.version)
    PythonCall.pyconvert_return(t(name, version))
end

function jl_convert_sim_circuit(t::Type{Braket.Circuit}, x)
    instructions = map(ix->pyconvert(Instruction, ix), x.instructions)
    results      = map(rt->pyconvert(Result, rt), x.results)
    prog         = t()
    foreach(ix->Braket.add_instruction!(prog, ix), instructions)
    foreach(rt->push!(prog.result_types, rt), results)
    PythonCall.pyconvert_return(prog)
end

py_obs(o::String) = pylist([pystr(o)])
function py_obs(obs::Vector)
    raw_obs = map(obs) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    return pylist(raw_obs)
end
Py(r::Braket.IR.Sample)        = braket[].ir.jaqcd.results.Sample(targets=convert_targets(r.targets), observable=py_obs(r.observable), type=pystr("sample"))
Py(r::Braket.IR.Expectation)   = braket[].ir.jaqcd.results.Expectation(targets=convert_targets(r.targets), observable=py_obs(r.observable), type=pystr("expectation"))
Py(r::Braket.IR.Variance)      = braket[].ir.jaqcd.results.Variance(observable=py_obs(r.observable), targets=convert_targets(r.targets), type=pystr("variance"))
Py(r::Braket.IR.Amplitude)     = braket[].ir.jaqcd.results.Amplitude(states=pylist(pystr(s) for s in r.states), type=pystr("amplitude"))
Py(r::Braket.IR.StateVector)   = braket[].ir.jaqcd.results.StateVector(type=pystr("statevector"))
Py(r::Braket.IR.DensityMatrix) = braket[].ir.jaqcd.results.DensityMatrix(targets=convert_targets(r.targets), type=pystr("densitymatrix"))
Py(r::Braket.IR.Probability)   = braket[].ir.jaqcd.results.Probability(targets=convert_targets(r.targets), type=pystr("probability"))
# exclude adjoint gradient translation from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
Py(r::Braket.IR.AdjointGradient) =  braket[].ir.jaqcd.results.AdjointGradient(targets=convert_targets(r.targets), observable=pylist(r.observable), parameters=pylist(pystr(p) for p in r.parameters), type=pystr("adjoint_gradient"))
# COV_EXCL_STOP

function Py(rt::Braket.ResultTypeValue)
    py_typ = Py(rt.type)
    py_val = if rt.value isa Dict
        pydict(rt.value)
    elseif rt.value isa Float64
        rt.value
    elseif rt.value isa Vector{Vector{Vector{Float64}}}
        pylist(pylist(pycomplex(v_...) for v_ in v) for v in rt.value)
    else
        pylist(rt.value)
    end
    return braket[].task_result.gate_model_task_result_v1.ResultTypeValue(type=py_typ, value=py_val)
end
function Py(ir::OpenQasmProgram)
    py_inputs = isnothing(ir.inputs) ? PythonCall.pybuiltins.None : pydict(ir.inputs)
    return braket[].ir.openqasm.program_v1.Program(source=pystr(ir.source), inputs=py_inputs)
end

function Py(r::GateModelTaskResult)
    py_measurements  = (isnothing(r.measurements) || isempty(r.measurements)) ? PythonCall.pybuiltins.None : pylist(pylist(meas) for meas in r.measurements)
    py_probabilities = !isnothing(r.measurementProbabilities) ? pydict(Dict(pystr(k)=>v for (k,v) in r.measurementProbabilities)) : PythonCall.pybuiltins.None
    py_qubits        = !isnothing(r.measuredQubits) ? pylist(r.measuredQubits) : PythonCall.pybuiltins.None
    py_results       = pylist(Py(rtv) for rtv in r.resultTypes)
    py_task_mtd      = braket[].task_result.task_metadata_v1.TaskMetadata(id=pystr(r.taskMetadata.id), shots=Py(r.taskMetadata.shots), deviceId=pystr(r.taskMetadata.deviceId))
    # we create a placeholder OpenQASM3 Program here -- this is replaced by the "true" action in the Python wrapper
    dummy_action     = braket[].ir.openqasm.program_v1.Program(source=pystr(""))
    py_addl_mtd      = braket[].task_result.additional_metadata.AdditionalMetadata(action=dummy_action)
    return braket[].task_result.GateModelTaskResult(measurements=py_measurements, measurementProbabilities=py_probabilities, resultTypes=py_results, measuredQubits=py_qubits, taskMetadata=py_task_mtd, additionalMetadata=py_addl_mtd)
end
