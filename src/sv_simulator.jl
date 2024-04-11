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
    ::Type{StateVectorSimulator{T,S}},
    qubit_count::Int,
) where {T,S<:AbstractVector{T}}
    sv = S(undef, 2^qubit_count)
    fill!(sv, zero(T))
    sv[1] = one(T)
    return sv
end

function StateVectorSimulator{T,S}(
    qubit_count::Int,
    shots::Int,
) where {T,S<:AbstractVector{T}}
    sv = init(StateVectorSimulator{T,S}, qubit_count)
    return StateVectorSimulator{T,S}(sv, qubit_count, shots)
end
"""
    StateVectorSimulator([::T], qubit_count::Int, shots::Int) -> StateVectorSimulator{T, Vector{T}}

Create a `StateVectorSimulator` with `2^qubit_count` elements and `shots` shots to be measured. The default element type is `ComplexF64`.
"""
StateVectorSimulator(::Type{T}, qubit_count::Int, shots::Int) where {T<:Number} =
    StateVectorSimulator{T,StateVector{T}}(qubit_count, shots)
StateVectorSimulator(qubit_count::Int, shots::Int) =
    StateVectorSimulator(ComplexF64, qubit_count, shots)
Braket.qubit_count(svs::StateVectorSimulator) = svs.qubit_count
"""
    properties(svs::StateVectorSimulator) -> GateModelSimulatorDeviceCapabilities

Query the properties and capabilities of a `StateVectorSimulator`, including which gates and result types are supported and the minimum and maximum shot and qubit counts.
"""
Braket.properties(svs::StateVectorSimulator) = sv_props
supported_operations(svs::StateVectorSimulator) =
    sv_props.action["braket.ir.openqasm.program"].supportedOperations
supported_result_types(svs::StateVectorSimulator) =
    sv_props.action["braket.ir.openqasm.program"].supportedResultTypes
Braket.device_id(svs::StateVectorSimulator) = "braket_sv_v2"
Braket.name(svs::StateVectorSimulator) = "StateVectorSimulator"
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

state_vector(svs::StateVectorSimulator) = svs.state_vector
density_matrix(svs::StateVectorSimulator) =
    kron(svs.state_vector, adjoint(svs.state_vector))

for (gate, obs) in (
    (:X, :(Braket.Observables.X)),
    (:Y, :(Braket.Observables.Y)),
    (:Z, :(Braket.Observables.Z)),
    (:I, :(Braket.Observables.I)),
    (:H, :(Braket.Observables.H)),
)
    @eval begin
        function apply_observable!(
            observable::$obs,
            sv::S,
            target::Int,
        ) where {S<:AbstractVector{<:Complex}}
            apply_gate!($gate(), sv, target)
            return sv
        end
    end
end
function apply_observable!(
    observable::Braket.Observables.HermitianObservable,
    sv::StateVector{T},
    target::Int,
) where {T<:Complex}
    n_amps = length(sv)
    mat = observable.matrix
    nq = Int(log2(n_amps))
    endian_qubit = nq - target - 1
    Threads.@threads for ix = 0:div(n_amps, 2)-1
        lower_ix = pad_bit(ix, endian_qubit)
        higher_ix = flip_bit(lower_ix, endian_qubit)
        ix_pair = [lower_ix + 1, higher_ix + 1]
        @views begin
            sv[ix_pair] = mat * sv[ix_pair]
        end
    end
    return sv
end
function apply_observable!(
    observable::Braket.Observables.HermitianObservable,
    sv::StateVector{T},
    t1::Int,
    t2::Int,
) where {T<:Complex}
    n_amps = length(sv)
    nq = Int(log2(n_amps))
    endian_t1 = nq - 1 - t1
    endian_t2 = nq - 1 - t2
    mat = observable.matrix
    small_t, big_t = minmax(endian_t1, endian_t2)
    Threads.@threads for ix = 0:div(n_amps, 4)-1
        # bit shift to get indices
        ix_00 = pad_bit(pad_bit(ix, small_t), big_t)
        ix_10 = flip_bit(ix_00, endian_t2)
        ix_01 = flip_bit(ix_00, endian_t1)
        ix_11 = flip_bit(ix_10, endian_t1)
        @views begin
            ind_vec = SVector{4,Int}(ix_00 + 1, ix_10 + 1, ix_01 + 1, ix_11 + 1)
            sv[ind_vec] = mat * sv[ind_vec]
        end
    end
    return sv
end

function apply_observable!(
    observable::Braket.Observables.HermitianObservable,
    sv::StateVector{T},
    t1::Int,
    t2::Int,
    targets::Int...,
) where {T<:Complex}
    n_amps = length(sv)
    nq = Int(log2(n_amps))
    ts = [t1, t2, targets...]
    endian_ts = nq - 1 .- ts
    o_mat = observable.matrix

    ordered_ts = sort(collect(endian_ts))
    flip_list = map(0:2^length(ts)-1) do t
        f_vals = Bool[(((1 << f_ix) & t) >> f_ix) for f_ix = 0:length(ts)-1]
        return ordered_ts[f_vals]
    end
    Threads.@threads for ix = 0:div(n_amps, 2^length(ts))-1
        padded_ix = ix
        for t in ordered_ts
            padded_ix = pad_bit(padded_ix, t)
        end
        ixs = map(flip_list) do f
            flipped_ix = padded_ix
            for f_val in f
                flipped_ix = flip_bit(flipped_ix, f_val)
            end
            return flipped_ix + 1
        end
        @views begin
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
    ev = dot(bra, ket)
    return real(ev)
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
    operations = reduce(vcat, diagonalizing_gates(obs...) for obs in observables)
    for op in operations
        apply_gate!(op.operator, svs._state_vector_after_observables, op.target...)
    end
    return svs
end
"""
    probabilities(svs::StateVectorSimulator) -> Vector{Float64}

Compute the observation probabilities of all amplitudes in the state vector in `svs`.
"""
probabilities(svs::StateVectorSimulator) = abs2.(svs.state_vector)

function _apply_ag_Hamiltonian(sv::StateVector{T}, H::Sum, targets) where {T<:Complex}
    bra = similar(sv)
    temp_sv = similar(sv)
    for (op, op_target) in zip(H.summands, targets)
        # don't conj here because H⁺ == H and sv will
        # be conj'd in dot
        copyto!(temp_sv, sv)
        apply_observable!(op, temp_sv, targets...)
        Threads.@threads for ix = 1:length(sv)
            bra[ix] += temp_sv[ix]
        end
    end
    return bra
end

function _apply_ag_Hamiltonian(
    sv::StateVector{T},
    H::Observable,
    targets,
) where {T<:Complex}
    bra = similar(sv)
    temp_sv = similar(sv)
    # don't conj here because H⁺ == H and sv will
    # be conj'd in dot
    copyto!(temp_sv, sv)
    apply_observable!(H, temp_sv, targets[1]...)
    Threads.@threads for ix = 1:length(sv)
        bra[ix] += temp_sv[ix]
    end
    return bra
end

function _get_params_and_angles(g::AngledGate{N}) where {N}
    params_and_angles = Dict{String,Int}()
    for (a_ix, angle) in enumerate(g.angle)
        angle isa FreeParameter && (params_and_angles[String(angle.name)] = a_ix)
    end
    return params_and_angles
end
_get_params_and_angles(g::Gate) = Dict{String,Int}()

function calculate(
    ag::AdjointGradient,
    instructions::Vector{Instruction},
    inputs::Dict{String,Float64},
    svs::StateVectorSimulator{T,S},
) where {T<:Complex,S<:AbstractStateVector{T}}
    svs.shots != 0 && throw(
        ArgumentError(
            "state vector simulator initiated with non-zero shots $(svs.shots). Adjoint gradient only supported for shots=0.",
        ),
    )
    # apply Hamiltonian to ket_svs
    H = ag.observable
    targets = ag.targets
    bra_sv = _apply_ag_Hamiltonian(svs.state_vector, H, targets)
    ket_sv = deepcopy(svs.state_vector)
    initial_ev = dot(bra_sv, ket_sv)
    rev_insts = Iterators.reverse(instructions)
    params_list = collect(keys(inputs))
    param_derivs =
        Dict{String,Float64}(zip(params_list, zeros(Float64, length(params_list))))
    # try to reuse this and avoid thrashing GC
    temp_ket_sv = deepcopy(svs.state_vector)
    for (ix, inst) in enumerate(rev_insts)
        params_and_angles = _get_params_and_angles(inst.operator)
        param_values =
            Dict{Symbol,Number}(Symbol(k) => inputs[k] for k in keys(params_and_angles))
        inv_gate = inverted_gate(bind_value!(inst.operator, param_values))
        apply_gate!(inv_gate, ket_sv, inst.target...)
        for (param, angle) in params_and_angles
            copyto!(temp_ket_sv, ket_sv)
            deriv_coeff, deriv_gate_fn = derivative_gate(
                Val(angle),
                bind_value!(inst.operator, param_values),
                inst.target...,
            )
            p_ix = findfirst(p == param for p in params_list)
            temp_ket_sv = deriv_gate_fn(temp_ket_sv)

            dot_prod = 2.0 * real(deriv_coeff * dot(bra_sv, temp_ket_sv))
            param_derivs[param] += dot_prod
        end
        if ix < length(instructions)
            # be careful about conjs here...
            apply_gate!(inv_gate, bra_sv, inst.target...)
        end
    end
    return initial_ev, param_derivs
end
