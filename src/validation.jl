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
    results,
    supported_result_types,
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
_validate_ir_results_compatibility(simulator::D, results, ::Val{V}) where {D<:AbstractSimulator, V} = _validate_ir_results_compatibility(simulator, results, supported_result_types(simulator, Val(V)))

function _validate_shots_and_ir_results(shots::Int, results, qubit_count::Int)
    if shots == 0
        (isnothing(results) || isempty(results)) && throw(ValidationError("Result types must be specified in the IR when shots=0", "ValueError"))
        for rt in results
            rt.type == "sample" && throw(ValidationError("sample can only be specified when shots>0", "ValueError"))
            rt.type == "amplitude" && _validate_amplitude_states(rt.states, qubit_count)
        end
    elseif shots > 0 && !isnothing(results) && !isempty(results)
        for rt in results
            rt.type ∈ ["statevector", "amplitude", "densitymatrix"] && throw(ValidationError(
                "statevector, amplitude and densitymatrix result types not available when shots>0",
                "ValueError"
            ))
        end
    end
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

function _validate_ir_instructions_compatibility(
    circuit::Union{Program,Circuit},
    supported_operations,
)
    circuit_instruction_names = map(ix->replace(lowercase(string(typeof(ix.operator))), "_"=>"", "braket."=>""), circuit.instructions)
    supported_instructions    = Set(map(op->replace(lowercase(op), "_"=>""), supported_operations))
    no_noise = true
    for name in circuit_instruction_names
        if name in _NOISE_INSTRUCTIONS
            no_noise = false
            if name ∉ supported_instructions
                throw(ValidationError("Noise instructions are not supported by the state vector simulator (by default). You need to use the density matrix simulator: LocalSimulator(\"braket_dm_v2\").", "TypeError"))
            end
        end
    end
    if no_noise && !isempty(intersect(_NOISE_INSTRUCTIONS, supported_instructions))
        @warn "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience."
    end
    return
end
_validate_ir_instructions_compatibility(simulator::D, circuit::Union{Program,Circuit}, v::Val{V}) where {D<:AbstractSimulator, V} = _validate_ir_instructions_compatibility(circuit, supported_operations(simulator, v))

_validate_result_type_qubits_exist(rt::Braket.StateVector, qubit_count::Int) = return
_validate_result_type_qubits_exist(rt::Braket.Amplitude, qubit_count::Int) = return
function _validate_result_type_qubits_exist(rt::AdjointGradient, qubit_count::Int)
    isempty(rt.targets) && return
    targets = reduce(vcat, targets)
    maximum(targets) > qubit_count &&
            throw(
                "Result type $rt references invalid qubits $targets. Maximum qubit number is $(qubit_count-1).",
            )
    return
end
# don't need to check for `isnothing` here as the `Braket.QubitSet` being empty covers this
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
    max_qc = 0
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
