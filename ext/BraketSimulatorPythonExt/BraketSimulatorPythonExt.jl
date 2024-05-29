module BraketSimulatorPythonExt

using BraketSimulator, BraketSimulator.Braket, PythonCall, BraketSimulator.Dates

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
    PythonCall.pyconvert_add_rule("braket.schema_common.schema_header:BraketSchemaHeader", Braket.braketSchemaHeader, jl_convert_bsh)
    PythonCall.pyconvert_add_rule("braket.ir.openqasm.program_v1:Program", OpenQasmProgram, jl_convert_oqprogram)
end

function simulate(
    simulator::AbstractSimulator,
    task_specs::Union{PyList{Any},NTuple{N,PyIterable}},
    args...;
    input::Union{PyList{Any},PyDict{Any,Any},Py}=PyDict{Any,Any}(),
    kwargs...,
) where {N}
    # handle inputs
    jl_inputs = nothing
    shots = args[end]
    stats = @timed begin
        if input isa PyDict{Any,Any}
            jl_inputs = pyconvert(Dict{String,Float64}, input)
        else
            jl_inputs = [pyconvert(Dict{String,Float64}, py_inputs) for py_inputs in input]
        end
        jl_specs = map(task_specs) do spec 
            return if pyhasattr(spec, "source")
                pyconvert(OpenQasmProgram, spec)
            else
                Braket.parse_raw_schema(pyconvert(String, spec.json()))
            end
        end
        input = nothing
    end
    @debug "Time for conversion of specs and inputs: $(stats.time)."
    if haskey(kwargs, :measured_qubits)
        jl_measured_qubits = [pyconvert(Int, qubit) for qubit in kwargs[:measured_qubits]]
        kwargs = merge(Dict(kwargs...), Dict(:measured_qubits=>jl_measured_qubits))
    end
    PythonCall.GC.disable()
    if length(jl_specs) == 1
        result     = simulate(simulator, jl_specs[1], args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        py_result  = Py(result)
    else # this is a batch! use a Braket.jl LocalSimulator to take advantage of thread migration
        local_sim   = Braket.LocalSimulator(simulator) 
        task_batch  = simulate(local_sim, jl_specs, args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        raw_results = results(task_batch)
        # now have to convert back to GateModelTaskResult from GateModelQuantumTaskResult
        processed_results = map(zip(raw_results, task_specs)) do (result, task_spec)
            header = Braket.braketSchemaHeader("braket.task_result.gate_model_task_result", "1")
            return Py(Braket.GateModelTaskResult(header, result.measurements, result.measurement_probabilities, result.result_types, result.measured_qubits, result.task_metadata, result.additional_metadata))
        end
        py_result = pylist(processed_results)
    end
    PythonCall.GC.enable()
    return py_result
end

end
