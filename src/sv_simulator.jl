"""
    StateVectorSimulator{T, S<:AbstractVector{T}} <: AbstractSimulator

Simulator representing a pure state evolution of a statevector of type `S`, with element type `T`. State vector simulators should be used to simulate circuits without noise.
"""
mutable struct StateVectorSimulator{T,S<:AbstractVector{T}} <: AbstractSimulator
    state_vector::S
    qubit_count::Int
    shots::Int
    buffer::AbstractArray{UInt8}
    shot_buffer::Vector{Int} # only used when shots>0
    _alias::Vector{Int}
    _ap::Vector{Float64}
    _larges::Vector{Int}
    _smalls::Vector{Int}
    _state_vector_after_observables::S
    function StateVectorSimulator{T,S}(
        state_vector::S,
        qubit_count::Int,
        shots::Int,
    ) where {T,S<:AbstractVector{T}}
        # careful here, need to dispatch on qubit count if it is large
        shot_buffer = Vector{Int}(undef, shots)
        ap_len  = ap_size(shots, qubit_count)
        _ap     = zeros(Float64, ap_len)
        _alias  = zeros(Int, ap_len)
        _larges = zeros(Int, ap_len)
        _smalls = zeros(Int, ap_len)
        return new(
            state_vector,
            qubit_count,
            shots,
            T[],
            shot_buffer,
            _alias,
            _ap,
            _larges,
            _smalls,
            S(undef, 0),
        )
    end
end
function init(
    t::Type{S},
    qubit_count::Int,
) where {T,S<:AbstractVector{T}}
    sv = t(undef, 2^qubit_count)
    fill!(sv, zero(T))
    sv[1] = one(T)
    return sv
end

function StateVectorSimulator{T,S}(
    qubit_count::Int,
    shots::Int,
) where {T,S<:AbstractVector{T}}
    sv = init(S, qubit_count)
    return StateVectorSimulator{T,S}(sv, qubit_count, shots)
end
"""
    StateVectorSimulator([::T], qubit_count::Int, shots::Int) -> StateVectorSimulator{T, Vector{T}}

Create a `StateVectorSimulator` with `2^qubit_count` elements and `shots` shots to be measured. The default element type is `ComplexF64`.
"""
StateVectorSimulator(::Type{T}, qubit_count::Int, shots::Int) where {T<:Number} =
    StateVectorSimulator{T,Vector{T}}(qubit_count, shots)
StateVectorSimulator(qubit_count::Int, shots::Int) =
    StateVectorSimulator(ComplexF64, qubit_count, shots)
qubit_count(svs::StateVectorSimulator) = svs.qubit_count
"""
    properties(svs::StateVectorSimulator) -> GateModelSimulatorDeviceCapabilities

Query the properties and capabilities of a `StateVectorSimulator`, including which gates and result types are supported and the minimum and maximum shot and qubit counts.
"""
properties(svs::StateVectorSimulator) = sv_props
supported_operations(svs::StateVectorSimulator, ::Val{:OpenQASM}) = sv_props.action["braket.ir.openqasm.program"].supportedOperations
supported_operations(svs::StateVectorSimulator) = supported_operations(svs::StateVectorSimulator, Val(:OpenQASM))
supported_operations(svs::StateVectorSimulator, ::Val{:JAQCD}) = sv_props.action["braket.ir.jaqcd.program"].supportedOperations
supported_result_types(svs::StateVectorSimulator, ::Val{:OpenQASM}) = sv_props.action["braket.ir.openqasm.program"].supportedResultTypes
supported_result_types(svs::StateVectorSimulator, ::Val{:JAQCD}) = sv_props.action["braket.ir.jaqcd.program"].supportedResultTypes
supported_result_types(svs::StateVectorSimulator) = supported_result_types(svs::StateVectorSimulator, Val(:OpenQASM))
device_id(svs::StateVectorSimulator) = "braket_sv_v2"
name(svs::StateVectorSimulator) = "StateVectorSimulator"
Base.show(io::IO, svs::StateVectorSimulator) =
    print(io, "StateVectorSimulator(qubit_count=$(qubit_count(svs)), shots=$(svs.shots))")
Base.similar(svs::StateVectorSimulator{T,S}; shots::Int = svs.shots) where {T,S} =
    StateVectorSimulator{T,S}(svs.qubit_count, shots)
Base.copy(svs::StateVectorSimulator{T,S}) where {T,S} =
    StateVectorSimulator{T,S}(deepcopy(svs.state_vector), svs.qubit_count, svs.shots)
function Base.copyto!(
    dst::StateVectorSimulator{T,S},
    src::StateVectorSimulator{T,S},
) where {T,S}
    copyto!(dst.state_vector, src.state_vector)
    return dst
end

function reinit!(
    svs::StateVectorSimulator{T,S},
    qubit_count::Int,
    shots::Int,
) where {T,S<:AbstractVector{T}}
    n = 2^qubit_count
    ap_len = ap_size(shots, qubit_count)
    if length(svs.state_vector) != n
        resize!(svs.state_vector, n)
        resize!(svs._alias,  ap_len)
        resize!(svs._ap,     ap_len)
        resize!(svs._larges, ap_len)
        resize!(svs._smalls, ap_len)
    end
    if svs.shots != shots
        resize!(svs._alias,  ap_len)
        resize!(svs._ap,     ap_len)
        resize!(svs._larges, ap_len)
        resize!(svs._smalls, ap_len)
        resize!(svs.shot_buffer, shots)
    end
    svs.state_vector .= zero(T)
    svs._ap          .= zero(Float64)
    svs._alias       .= zero(Int)
    svs._larges      .= zero(Int)
    svs._smalls      .= zero(Int)
    svs.state_vector[1] = one(T)
    svs.qubit_count = qubit_count
    svs.shots = shots
    svs._state_vector_after_observables = S(undef, 0)
    return
end

"""
    evolve!(svs::StateVectorSimulator{T, S<:AbstractVector{T}}, operations::Vector{Instruction}) -> StateVectorSimulator{T, S}

Apply each operation of `operations` in-place to the state vector contained in `svs`.

Effectively, perform the operation:

`` \\left| \\psi \\right\\rangle \\to \\hat{A} \\left| \\psi \\right\\rangle ``

for each operation ``\\hat{A}`` in `operations`.
"""
function evolve!(
    svs::StateVectorSimulator{T,S},
    operations,
)::StateVectorSimulator{T,S} where {T<:Complex,S<:AbstractVector{T}}
    for (oix, op) in enumerate(operations)
        apply_gate!(op.operator, svs.state_vector, op.target...)
    end
    return svs
end

state_vector(svs::StateVectorSimulator) = copy(svs.state_vector)
density_matrix(svs::StateVectorSimulator) =
    kron(svs.state_vector, adjoint(svs.state_vector))

function apply_observable!(
    gate::G,
    sv::S,
    target::Int,
) where {S<:AbstractVector{<:Complex}, G<:Gate}
    apply_gate!(gate, sv, target)
    return sv
end
function apply_observable!(
    observable::Observables.HermitianObservable,
    sv::AbstractStateVector{T},
    target::Int,
) where {T<:Complex}
    n_amps   = length(sv)
    mat      = observable.matrix
    n_qubits = Int(log2(n_amps))
    endian_qubit = n_qubits - target - 1
    Threads.@threads for ix = 0:div(n_amps, 2)-1
        lower_ix  = pad_bit(ix, endian_qubit)
        higher_ix = flip_bit(lower_ix, endian_qubit)
        ix_pair   = [lower_ix + 1, higher_ix + 1]
        @views @inbounds sv[ix_pair] = mat * sv[ix_pair]
    end
    return sv
end
function apply_observable!(
    observable::Observables.HermitianObservable,
    sv::AbstractStateVector{T},
    target_1::Int,
    target_2::Int,
) where {T<:Complex}
    n_amps, (endian_t1, endian_t2) = get_amps_and_qubits(sv, target_1, target_2)
    small_t, big_t = minmax(endian_t1, endian_t2)
    mat            = observable.matrix
    Threads.@threads for ix = 0:div(n_amps, 4)-1
        # bit shift to get indices
        ix_00 = pad_bits(ix, (small_t, big_t))
        ix_10 = flip_bit(ix_00, endian_t2)
        ix_01 = flip_bit(ix_00, endian_t1)
        ix_11 = flip_bit(ix_10, endian_t1)
        ind_vec = SVector{4,Int}(ix_00 + 1, ix_10 + 1, ix_01 + 1, ix_11 + 1)
        @views @inbounds sv[ind_vec] = mat * sv[ind_vec]
    end
    return sv
end

function apply_observable!(
    observable::Observables.HermitianObservable,
    sv::AbstractStateVector{T},
    target_1::Int,
    target_2::Int,
    targets::Int...,
) where {T<:Complex}
    n_amps   = length(sv)
    n_qubits = Int(log2(n_amps))
    full_targets = [target_1, target_2, targets...]
    endian_ts = n_qubits - 1 .- full_targets
    o_mat     = observable.matrix
    n_targets = length(full_targets)

    ordered_ts = sort(collect(endian_ts))
    flip_list  = map(0:2^n_targets-1) do target_amp
        to_flip = Bool[(((1 << target_ix) & target_amp) >> target_ix) for target_ix = 0:n_targets-1]
        return ordered_ts[to_flip]
    end
    Threads.@threads for ix = 0:div(n_amps, 2^n_targets)-1
        padded_ix = pad_bits(ix, ordered_ts)
        ixs = map(flip_list) do bits_to_flip
            flip_bits(padded_ix, bits_to_flip) + 1
        end
        @views @inbounds begin
            amps = sv[ixs[:]]
            sv[ixs[:]] = o_mat * amps
        end
    end
    return sv
end

"""
    expectation(svs::StateVectorSimulator, observable::Observables.Observable, targets::Int...) -> Float64

Compute the exact (`shots=0`) expectation value of `observable` applied to `targets`
given the evolved state vector in `svs`. In other words, compute

``\\langle \\psi | \\hat{O} | \\psi \\rangle``.
"""
function expectation(
    svs::StateVectorSimulator,
    observable::Observables.Observable,
    targets::Int...,
)
    bra = svs.state_vector
    ket = apply_observable(observable, svs.state_vector, targets...)
    expectation_value = dot(bra, ket)
    return real(expectation_value)
end

function state_with_observables(svs::StateVectorSimulator)
    isempty(svs._state_vector_after_observables) &&
        error("observables have not been applied.")
    return svs._state_vector_after_observables
end

function apply_observables!(svs::StateVectorSimulator, observables)
    !isempty(svs._state_vector_after_observables) &&
        error("observables have already been applied.")
    svs._state_vector_after_observables = deepcopy(svs.state_vector)
    operations = mapreduce(obs->diagonalizing_gates(obs...), vcat, observables)
    for operation in operations
        apply_gate!(operation.operator, svs._state_vector_after_observables, operation.target...)
    end
    return svs
end
"""
    probabilities(svs::StateVectorSimulator) -> Vector{Float64}

Compute the observation probabilities of all amplitudes in the state vector in `svs`.
"""
probabilities(svs::StateVectorSimulator) = abs2.(svs.state_vector)
