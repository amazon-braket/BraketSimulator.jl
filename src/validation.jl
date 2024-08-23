struct ValidationError <: Exception
    message::String
    alternate_type::String
end

const _NOISE_INSTRUCTIONS = Set(replace(lowercase(op), "_"=>"") for op in
                                ["amplitude_damping",
                                 "bit_flip",
                                 "depolarizing",
                                 "generalized_amplitude_damping",
                                 "kraus",
                                 "pauli_channel",
                                 "two_qubit_pauli_channel",
                                 "phase_flip",
                                 "phase_damping",
                                 "two_qubit_dephasing",
                                 "two_qubit_depolarizing",
                                ]
                            )

function _validate_amplitude_states(states::Vector{String}, qubit_count::Int)
    !all(length(state) == qubit_count for state in states) && throw(ValidationError("all states in $states must have length $qubit_count.", "ValueError"))
    return
end

function _validate_ir_results_compatibility(
    simulator::D,
    results::Vector{AbstractProgramResult},
    supported_result_types::Vector,
) where {D<:AbstractSimulator}
    (isnothing(results) || isempty(results)) && return
    circuit_result_types_name    = [rt.type for rt in results]
    supported_result_types_names = map(supported_rt->lowercase(string(supported_rt.name)), supported_result_types)
    for name in circuit_result_types_name
        name ∉ supported_result_types_names &&
            throw(ErrorException("result type $name not supported by $simulator"))
    end
    return
end
function _validate_ir_results_compatibility(
    simulator::D,
    results::Vector{Result},
    supported_result_types::Vector,
) where {D<:AbstractSimulator}
    (isnothing(results) || isempty(results)) && return
    circuit_result_types_name    = map(label, results)
    supported_result_types_names = map(supported_rt->lowercase(string(supported_rt.name)), supported_result_types)
    for name in circuit_result_types_name
        name ∉ supported_result_types_names &&
            throw(ErrorException("result type $name not supported by $simulator"))
    end
    return
end
_validate_ir_results_compatibility(simulator::D, results::Vector{A}, ::Val{V}) where {D<:AbstractSimulator, A, V} = _validate_ir_results_compatibility(simulator, results, supported_result_types(simulator, Val(V)))

_validate_exact_result(amp::Amplitude, qubit_count)    = _validate_amplitude_states(amp.states, qubit_count)
_validate_exact_result(amp::IR.Amplitude, qubit_count) = _validate_amplitude_states(amp.states, qubit_count)
_validate_exact_result(sample::Sample, qubit_count)    = throw(ValidationError("sample can only be specified when shots>0", "ValueError"))
_validate_exact_result(sample::IR.Sample, qubit_count) = throw(ValidationError("sample can only be specified when shots>0", "ValueError"))
_validate_exact_result(result, qubit_count)            = return

const nonzero_shots_error = "statevector, amplitude and densitymatrix result types are only available when shots==0"
_validate_nonzero_shots_result(sv::StateVector)      = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(sv::IR.StateVector)   = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(amp::Amplitude)       = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(amp::IR.Amplitude)    = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(dm::DensityMatrix)    = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(dm::IR.DensityMatrix) = throw(ValidationError(nonzero_shots_error, "ValueError"))
_validate_nonzero_shots_result(result)               = return
function _validate_shots_and_ir_results(shots::Int, results, qubit_count::Int)
    if shots == 0
        (isnothing(results) || isempty(results)) && throw(ValidationError("Result types must be specified in the IR when shots=0", "ValueError"))
        foreach(r->_validate_exact_result(r, qubit_count), results)
    elseif shots > 0 && !isnothing(results)
        foreach(_validate_nonzero_shots_result, results)
    end
    return
end

function _validate_input_provided(circuit)
    for instruction in circuit.instructions
        possible_parameters = (Symbol("_angle"), Symbol("_angle_1"), Symbol("_angle_2"), Symbol("angle"))
        for parameter_name in possible_parameters
            if hasproperty(instruction.operator, parameter_name)
                for param in getproperty(instruction.operator, parameter_name)
                    try
                        Float64(param)
                    catch ex
                        throw(ErrorException("Missing input variable '$param'."))
                    end
                end
            end
        end
    end
    return
end

function _validate_ir_instructions_compatibility(circuit, supported_operations)
    circuit_instruction_names = map(ix->replace(lowercase(string(typeof(ix.operator))), "_"=>"", "braketsimulator."=>""), circuit.instructions)
    supported_instructions    = Set(map(op->replace(lowercase(op), "_"=>""), supported_operations))
    no_noise = true
    for name in filter(ix_name->ix_name ∈ _NOISE_INSTRUCTIONS, circuit_instruction_names)
        no_noise = false
        name ∉ supported_instructions && throw(ValidationError("Noise instructions are not supported by the state vector simulator (by default). You need to use the density matrix simulator: LocalSimulator(\"braket_dm_v2\").", "TypeError"))
    end
    if no_noise && !isempty(intersect(_NOISE_INSTRUCTIONS, supported_instructions))
        @warn "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience."
    end
    return
end
_validate_ir_instructions_compatibility(simulator::D, circuit, v::Val{V}) where {D<:AbstractSimulator, V} = _validate_ir_instructions_compatibility(circuit, supported_operations(simulator, v))

_validate_result_type_qubits_exist(rt::StateVector, qubit_count::Int) = return
_validate_result_type_qubits_exist(rt::Amplitude, qubit_count::Int)   = return
# exclude adjoint gradient validation from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
function _validate_result_type_qubits_exist(rt::AdjointGradient, qubit_count::Int)
    isempty(rt.targets) && return
    targets = reduce(vcat, targets)
    maximum(targets) > qubit_count &&
            throw(
                "Result type $rt references invalid qubits $targets. Maximum qubit number is $(qubit_count-1).",
            )
    return
end
# COV_EXCL_STOP
# don't need to check for `isnothing` here as the `QubitSet` being empty covers this
function _validate_result_type_qubits_exist(rt::RT, qubit_count::Int) where {RT<:Result}
    isempty(rt.targets) && return
    maximum(rt.targets) > qubit_count &&
            throw(
                  "Result type $rt references invalid qubits $(rt.targets). Maximum qubit number is $(qubit_count-1).",
            )
    return
end
_validate_result_types_qubits_exist(rts::Vector{RT}, qubit_count::Int) where {RT<:Result} = foreach(rt->_validate_result_type_qubits_exist(rt, qubit_count), rts)

function _validate_operation_qubits(operations::Vector{<:Instruction})
    unique_qs = Set{Int}()
    max_qc    = 0
    for operation in operations
        targets = operation.target
        max_qc = max(max_qc, targets...)
        union!(unique_qs, targets)
    end
    max_qc >= length(unique_qs) && throw(ValidationError(
        "Non-contiguous qubit indices supplied; qubit indices in a circuit must be contiguous. Qubits referenced: $unique_qs",
        "ValueError"
    ))
    return
end

function _check_observable(observable_map, observable, qubits)
    for qubit in qubits, (target, previously_measured) in observable_map
        qubit ∉ target && continue
        # must ensure that target is the same
        issetequal(target, qubits) || throw(ValidationError("Qubit part of incompatible results targets", "ValueError"))
        # must ensure observable is the same
        typeof(previously_measured) != typeof(observable) && throw(ValidationError("Conflicting result types applied to a single qubit", "ValueError"))
        # including matrix value for Hermitians
        observable isa Observables.HermitianObservable && !isapprox(previously_measured.matrix, observable.matrix) && throw(ValidationError("Conflicting result types applied to a single qubit", "ValueError"))
    end
    observable_map[qubits] = observable
    return observable_map
end

function _combine_obs_and_targets(observable::Observables.HermitianObservable, result_targets::Vector{Int})
    obs_qc = qubit_count(observable)
    length(result_targets) == obs_qc && return [(observable, result_targets)]
    obs_qc == 1 && return [(copy(observable), rt) for rt in result_targets]
    throw(ValidationError("only single-qubit Hermitian observables can be applied to all qubits.", "ValueError"))
end
_combine_obs_and_targets(observable::Observables.TensorProduct, result_targets::Vector{Int}) = zip(observable.factors, [splice!(result_targets, 1:qubit_count(f)) for f in observable.factors])
_combine_obs_and_targets(observable, result_targets::Vector{Int}) = length(result_targets) == 1 ? [(observable, result_targets)] : [(copy(observable), t) for t in result_targets]

function _verify_openqasm_shots_observables(circuit::Circuit, n_qubits::Int)
    observable_map = Dict()
    for result in filter(rt->rt isa ObservableResult, circuit.result_types)
        result.observable isa Observables.I && continue
        result_targets = isempty(result.targets) ? collect(0:n_qubits-1) : collect(result.targets)
        for obs_and_target in _combine_obs_and_targets(result.observable, result_targets)
            observable_map = _check_observable(observable_map, obs_and_target...)
        end
    end
    return
end
