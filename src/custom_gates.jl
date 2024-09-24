mutable struct DoubleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    DoubleExcitation(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{DoubleExcitation}) = 4
function matrix_rep_raw(::DoubleExcitation, ϕ) # nosemgrep
    sθ, cθ      = sincos(ϕ/2.0)
    mat         = diagm(ones(ComplexF64, 16))
    mat[4, 4]   = cθ
    mat[13, 13] = cθ
    mat[4, 13]  = -sθ
    mat[13, 4]  = sθ
    return SMatrix{16,16,ComplexF64}(mat)
end

mutable struct SingleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    SingleExcitation(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{SingleExcitation}) = 2
matrix_rep_raw(::SingleExcitation, ϕ) = ((sθ, cθ) = sincos(ϕ/2.0); return SMatrix{4,4,ComplexF64}(complex(1.0), 0, 0, 0, 0, cθ, -sθ, 0, 0, sθ, cθ, 0, 0, 0, 0, complex(1.0)))
"""
    MultiRz(angle)

Multi-qubit Z-rotation gate. The 2-qubit version is equivalent to the `ZZ` gate,
and the single-qubit version is equivalent to the `Rz` gate.
"""
struct MultiRZ <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    MultiRZ(angle::T, pow_exponent::Float64=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} = new(angle, pow_exponent)
end
qubit_count(g::MultiRZ) = 1

"""
    MultiQubitPhaseShift{N}(angle)

Global phase shift on `N` qubits. Equivalent to
the OpenQASM3 built-in [`gphase` gate](https://openqasm.com/language/gates.html#gphase).
Controls/negative controls applied to this gate control
which states are rotated, so that `Control(g::MultiQubitPhaseShift{2})` will apply the rotation
to the `|11>` state.
"""
mutable struct MultiQubitPhaseShift{N} <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    MultiQubitPhaseShift{N}(angle::T, pow_exponent::Float64=1.0) where {N, T<:NTuple{1,Union{Real,FreeParameter}}} = new(angle, pow_exponent)
end
matrix_rep_raw(g::MultiQubitPhaseShift{N}) where {N} = Diagonal(SVector{2^N, ComplexF64}(exp(im*g.angle[1]) for _ in 1:2^N))
qubit_count(g::MultiQubitPhaseShift{N}) where {N} = N

function apply_gate!(
    factor::ComplexF64,
    diagonal::PauliEigenvalues{N},
    state_vec::AbstractStateVector{T},
    ts::Vararg{Int,N},
) where {T<:Complex,N}
    g_mat  = Diagonal(SVector{2^N,ComplexF64}(exp(factor * diagonal[i]) for i = 1:2^N))
    apply_gate!(g_mat, state_vec, ts...)
end
function apply_gate!(
    factor::ComplexF64,
    state_vec::AbstractStateVector{T},
    ts::Vararg{Int,N},
) where {T<:Complex,N}
    g_mat = Diagonal(SVector{2^N,ComplexF64}(factor for i = 1:2^N))
    apply_gate!(g_mat, state_vec, ts...)
end
for (V, f) in ((true, :conj), (false, :identity))
    @eval begin
        apply_gate!(::Val{$V}, g::MultiRZ, state_vec::AbstractStateVector{T}, t::Int) where {T<:Complex} = apply_gate!(Val($V), Rz(g.angle, g.pow_exponent), state_vec, t)
        apply_gate!(::Val{$V}, g::MultiRZ, state_vec::AbstractStateVector{T}, t1::Int,t2::Int,) where {T<:Complex} = apply_gate!(Val($V), ZZ(g.angle, g.pow_exponent), state_vec, t1, t2)
        apply_gate!(::Val{$V}, g::MultiRZ, state_vec::AbstractStateVector{T}, ts::Vararg{Int,N}) where {T<:Complex,N} = apply_gate!($f(-im * g.angle[1] / 2.0), PauliEigenvalues(Val(N)), state_vec, ts...)
        apply_gate!(::Val{$V}, g::MultiQubitPhaseShift{N}, state_vec::AbstractStateVector{T}, ts::Vararg{Int,N}) where {T<:Complex, N} = apply_gate!($f(exp(im*g.angle[1]*g.pow_exponent)), state_vec, ts...)
    end
end

function apply_gate!(
    g::DoubleExcitation,
    state_vec::AbstractStateVector{T},
    t1::Int,
    t2::Int,
    t3::Int,
    t4::Int,
) where {T<:Complex}
    n_amps, endian_ts = get_amps_and_qubits(state_vec, t1, t2, t3, t4)
    ordered_ts = sort(collect(endian_ts))
    sinϕ, cosϕ = sincos(g.angle[1] * g.pow_exponent / 2.0)
    e_t1, e_t2, e_t3, e_t4 = endian_ts
    Threads.@threads for ix = 0:div(n_amps, 2^4)-1
        padded_ix = pad_bits(ix, ordered_ts)
        i0011     = flip_bits(padded_ix, (e_t3, e_t4)) + 1
        i1100     = flip_bits(padded_ix, (e_t1, e_t2)) + 1
        @inbounds begin
            amp0011 = state_vec[i0011]
            amp1100 = state_vec[i1100]
            state_vec[i0011] = cosϕ * amp0011 - sinϕ * amp1100
            state_vec[i1100] = sinϕ * amp0011 + cosϕ * amp1100
        end
    end
    return
end
apply_gate!(::Val{true}, g::DoubleExcitation, state_vec::AbstractStateVector{T}, t1::Int, t2::Int, t3::Int, t4::Int) where {T<:Complex} = apply_gate!(g, state_vec, t1, t2, t3, t4)
apply_gate!(::Val{false}, g::DoubleExcitation, state_vec::AbstractStateVector{T}, t1::Int, t2::Int, t3::Int, t4::Int) where {T<:Complex} = apply_gate!(g, state_vec, t1, t2, t3, t4)

mutable struct Control{G<:Gate, B} <: Gate
    g::G
    bitvals::NTuple{B, Int}
    pow_exponent::Float64
    Control{G, B}(g::G, bitvals::NTuple{B, Int}, pow_exponent::Float64=1.0) where {B, G} = new(g, bitvals, pow_exponent)
end
Control(g::G, bitvals::NTuple{B, Int}, pow_exponent::Float64=1.0) where {G<:Gate, B} = Control{G, B}(g, bitvals, pow_exponent)
Control(g::Control{G, BC}, bitvals::NTuple{B, Int}, pow_exponent::Float64=1.0) where {G<:Gate, BC, B} = Control(g.g, (g.bitvals..., bitvals...), g.pow_exponent * pow_exponent)
Control(g::Control{G, BC}, bitvals::NTuple{0, Int}, pow_exponent::Float64=1.0) where {G<:Gate, BC} = g
Control(g::G, bitvals::NTuple{0, Int}, pow_exponent::Float64=1.0) where {G<:Gate} = g
Base.copy(c::Control{G, B}) where {G<:Gate, B} = Control(copy(c.g), c.bitvals, c.pow_exponent)
qubit_count(c::Control{G, B}) where {G<:Gate, B} = qubit_count(c.g) + B
qubit_count(c::Control{MultiQubitPhaseShift{N}, B}) where {N, B} = N
function matrix_rep_raw(c::Control{G, B}) where {G<:Gate, B}
    inner_mat = matrix_rep(c.g)
    inner_qc  = qubit_count(c.g)
    total_qc  = qubit_count(c.g) + B
    full_mat  = diagm(ones(ComplexF64, 2^total_qc))
    ctrl_ix   = 0
    for (b_ix, b) in enumerate(c.bitvals)
        ctrl_ix += b << (b_ix + inner_qc - 1)
    end
    for inner_ix in 1:2^qubit_count(c.g), inner_jx in 1:2^qubit_count(c.g)
        full_mat[ctrl_ix + inner_ix, ctrl_ix + inner_jx] = inner_mat[inner_ix, inner_jx]
    end
    return full_mat
end
function matrix_rep_raw(c::Control{MultiQubitPhaseShift{N}, B}) where {N, B}
    g_mat    = matrix_rep(c.g)
    qc       = N 
    ctrl_ix  = 0
    full_mat = ones(ComplexF64, 2^N)
    for (b_ix, b) in enumerate(c.bitvals)
        ctrl_ix += b << (b_ix - 1)
    end
    for inner_ix in 1:2^N
        full_mat[inner_ix] = (inner_ix == ctrl_ix+1) ? g_mat[inner_ix, inner_ix] : one(ComplexF64)
    end
    return Diagonal(SVector{2^N, ComplexF64}(full_mat))
end
