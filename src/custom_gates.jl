"""
    DoubleExcitation(ϕ)

Generate the matrix representation of the [DoubleExcitation](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitation.html) gate.

This gate performs an SO(2) rotation in the subspace ``{|1100\\rangle, |0011\\rangle}``, transforming the states as follows:

```math
|0011\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|0011\\rangle + \\sin\\left(\\frac{\\phi}{2}\\right)|1100\\rangle \\
|1100\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|1100\\rangle - \\sin\\left(\\frac{\\phi}{2}\\right)|0011\\rangle
```

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix =  BraketSimulator.DoubleExcitation(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

julia> eq2 = m * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];

julia> eq1 ==  [0, 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, 0] == true;

julia> eq2 ==  [0, 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, 0] == true;
```

"""
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

"""
    DoubleExcitationPlus(ϕ)

Generate the matrix representation of the [DoubleExcitationPlus](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitationPlus.html) gate.

This gate performs an SO(2) rotation in the subspace ``{|1100\\rangle, |0011\\rangle}`` with a phase-shift on other states:

```math
|0011\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|0011\\rangle - \\sin\\left(\\frac{\\phi}{2}\\right)|1100\\rangle \\
|1100\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|1100\\rangle + \\sin\\left(\\frac{\\phi}{2}\\right)|0011\\rangle \\
|x\\rangle & \\rightarrow e^{\\frac{i\\phi}{2}}|x\\rangle \\quad \\text{for all other basis states } |x\\rangle
```

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.DoubleExcitationPlus(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

julia> eq2 = m * [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];

julia> eq1 ==  [exp(im*ϕ/2), 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, exp(im*ϕ/2)] == true;

julia> eq2 ==  [exp(im*ϕ/2), 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, exp(im*ϕ/2)] == true;
```

"""
struct DoubleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    DoubleExcitationPlus(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{DoubleExcitationPlus}) = 4
function matrix_rep_raw(::DoubleExcitationPlus, ϕ) # nosemgrep
    sϕ, cϕ = sincos(ϕ / 2.0)
    eiϕ2   = exp(im * ϕ / 2.0)
    mat    = diagm(eiϕ2 * ones(ComplexF64, 16))
    @views begin
        mat[4, :]  .= 0
        mat[:, 4]  .= 0
        mat[13, :] .= 0
        mat[:, 13] .= 0
    end
    # Apply phase-shift to states outside rotation subspace
    mat[4, 4]   = cϕ
    mat[13, 13] = cϕ
    mat[4, 13]  = -sϕ
    mat[13, 4]  = sϕ
    return SMatrix{16, 16, ComplexF64}(mat)
end

"""
    DoubleExcitationMinus(ϕ)

Generate the matrix representation of the [DoubleExcitationMinus](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitationMinus.html) gate.

This gate performs an SO(2) rotation in the subspace ``{|1100\\rangle, |0011\\rangle}`` with a phase-shift on other states:

```math
|0011\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|0011\\rangle - \\sin\\left(\\frac{\\phi}{2}\\right)|1100\\rangle \\
|1100\\rangle & \\rightarrow \\cos\\left(\\frac{\\phi}{2}\\right)|1100\\rangle + \\sin\\left(\\frac{\\phi}{2}\\right)|0011\\rangle \\
|x\\rangle & \\rightarrow e^{-\\frac{i\\phi}{2}}|x\\rangle \\quad \\text{for all other basis states } |x\\rangle
```

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.DoubleExcitationMinus(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

julia> eq2 = m * [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];

julia> eq1 ==  [exp(-im*ϕ/2), 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, exp(-im*ϕ/2)] == true;

julia> eq2 ==  [exp(-im*ϕ/2), 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, exp(-im*ϕ/2)] == true;
```

"""
struct DoubleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    DoubleExcitationMinus(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{DoubleExcitationMinus}) = 4
function matrix_rep_raw(::DoubleExcitationMinus, ϕ) # nosemgrep
    sϕ, cϕ = sincos(ϕ / 2.0)
    eiϕ2   = exp(-im * ϕ / 2.0)
    mat    = diagm(eiϕ2 * ones(ComplexF64, 16))
    @views begin
        mat[4, :]  .= 0
        mat[:, 4]  .= 0
        mat[13, :] .= 0
        mat[:, 13] .= 0
    end
    # Apply phase-shift to states outside rotation subspace
    mat[4, 4]   = cϕ
    mat[13, 13] = cϕ
    mat[4, 13]  = -sϕ
    mat[13, 4]  = sϕ
    return SMatrix{16, 16, ComplexF64}(mat)
end

"""
    SingleExcitation(ϕ)

Generate the matrix representation of the [SingleExcitation](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitation.html) gate. 

This gate performs a rotation in the subspace ``{|01\\rangle, |10\\rangle}``.

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.SingleExcitation(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [0, 1, 0, 0];

julia> eq2 = m * [0, 0, 1, 0];

julia> eq1 == [0, cos(ϕ/2), - sin(ϕ/2), 0] == true;
 
julia> eq2 == [0, sin(ϕ/2), cos(ϕ/2), 0] == true;
```

"""
mutable struct SingleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    SingleExcitation(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{SingleExcitation}) = 2
matrix_rep_raw(::SingleExcitation, ϕ) = ((sθ, cθ) = sincos(ϕ/2.0); return SMatrix{4,4,ComplexF64}(complex(1.0), 0, 0, 0, 0, cθ, -sθ, 0, 0, sθ, cθ, 0, 0, 0, 0, complex(1.0)))
"""
    SingleExcitationPlus(ϕ)

Generate the matrix representation of the [SingleExcitationPlus](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitationPlus.html) gate.

This gate performs a rotation in the subspace ``{|01\\rangle, |10\\rangle}`` with a phase-shift.

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.SingleExcitationPlus(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [1, 1, 0, 0];

julia> eq2 = m * [1, 0, 1, 0];

julia> eq1 == [exp(im*ϕ/2), cos(ϕ/2), - sin(ϕ/2), 0] == true; 

julia> eq2 == [exp(im*ϕ/2), sin(ϕ/2), cos(ϕ/2), 0] == true;
```

"""
struct SingleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    SingleExcitationPlus(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{SingleExcitationPlus}) = 2
function matrix_rep_raw(::SingleExcitationPlus, ϕ) # nosemgrep
    sϕ, cϕ = sincos(ϕ / 2.0)
    eiϕ2   = exp(im * ϕ / 2.0)
    return SMatrix{4,4,ComplexF64}(eiϕ2, 0, 0, 0, 0, cϕ, -sϕ, 0, 0, sϕ, cϕ, 0, 0, 0, 0, eiϕ2)
end

"""
    SingleExcitationMinus(ϕ)

Generate the matrix representation of the [SingleExcitationMinus](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitationMinus.html) gate.

This gate performs a rotation in the subspace ``{|01\\rangle, |10\\rangle}`` with a phase-shift.

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.SingleExcitationMinus(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [1, 1, 0, 0];

julia> eq2 = m * [1, 0, 1, 0];

julia> eq1 == [exp(-im*ϕ/2), cos(ϕ/2), - sin(ϕ/2), 0] == true;

julia> eq2 == [exp(-im*ϕ/2), sin(ϕ/2), cos(ϕ/2), 0] == true;
```

"""
struct SingleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    SingleExcitationMinus(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{SingleExcitationMinus}) = 2
function matrix_rep_raw(::SingleExcitationMinus, ϕ) # nosemgrep
    sϕ, cϕ = sincos(ϕ / 2.0)
    eiϕ2   = exp(-im * ϕ / 2.0)
    return SMatrix{4,4,ComplexF64}(eiϕ2, 0, 0, 0, 0, cϕ, -sϕ, 0, 0, sϕ, cϕ, 0, 0, 0, 0, eiϕ2)
end

"""
    FermionicSWAP(ϕ)

Generate the matrix representation of the [FermionicSWAP](https://docs.pennylane.ai/en/stable/code/api/pennylane.FermionicSWAP.html) gate.

This gate performs a rotation in adjacent fermionic modes under the Jordan-Wigner mapping, transforming states as follows:

```math
|00\\rangle & \\rightarrow |00\\rangle \\
|01\\rangle & \\rightarrow e^{\\frac{i\\phi}{2}}\\cos\\left(\\frac{\\phi}{2}\\right)|01\\rangle - ie^{\\frac{i\\phi}{2}}\\sin\\left(\\frac{\\phi}{2}\\right)|10\\rangle \\
|10\\rangle & \\rightarrow -ie^{\\frac{i\\phi}{2}}\\sin\\left(\\frac{\\phi}{2}\\right)|01\\rangle + e^{\\frac{i\\phi}{2}}\\cos\\left(\\frac{\\phi}{2}\\right)|10\\rangle \\
|11\\rangle & \\rightarrow e^{i\\phi}|11\\rangle
```

# Examples

```jldoctest
julia> using BraketSimulator

julia> ϕ = 3.56;

julia> gate_matrix = BraketSimulator.FermionicSWAP(ϕ);

julia> m  = BraketSimulator.matrix_rep(gate_matrix);

julia> eq1 = m * [0, 0, 0, 0];

julia> eq2 = m * [0, 1, 0, 0];

julia> eq3 = m * [0, 0, 1, 0];

julia> eq4 = m * [0, 0, 0, 1];

julia> eq1 == [0, 0, 0, 0] == true;
 
julia> eq2 == [0, exp(im*ϕ/2.0)*cos(ϕ / 2.0), - im*exp(im*ϕ/2.0)*sin(ϕ/2.0), 0] == true;

julia> eq3 == [0, - im*exp(im*ϕ/2.0)*sin(ϕ/2.0), exp(im*ϕ/2.0)*cos(ϕ/2.0), 0] == true;

julia> eq4 == [0, 0, 0, exp(im * ϕ)] == true;
```

"""
struct FermionicSWAP <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    FermionicSWAP(angle::T, pow_exponent=1.0) where {T<:NTuple{1,Union{Real,FreeParameter}}} =
        new(angle, Float64(pow_exponent))
end
qubit_count(::Type{FermionicSWAP}) = 2
function matrix_rep_raw(::FermionicSWAP, ϕ) # nosemgrep
    sϕ, cϕ = sincos(ϕ / 2.0)
    eiϕ2   = exp(im * ϕ / 2.0)
    eiϕ    = exp(im * ϕ)
    ieiϕ2  = im * eiϕ2
    return SMatrix{4,4,ComplexF64}(1, 0, 0, 0, 0, eiϕ2 * cϕ, -ieiϕ2 * sϕ, 0, 0, -ieiϕ2*sϕ, eiϕ2*cϕ, 0, 0, 0, 0, eiϕ)
end

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
    GPhase{N}(angle)

Global phase shift on `N` qubits. Equivalent to
the OpenQASM3 built-in [`gphase` gate](https://openqasm.com/language/gates.html#gphase).
Controls/negative controls applied to this gate control
which states are rotated, so that `Control(g::GPhase{2})` will apply the rotation
to the `|11>` state.
"""
mutable struct GPhase{N} <: AngledGate{1}
    angle::NTuple{1,Union{Real,FreeParameter}}
    pow_exponent::Float64
    GPhase{N}(angle::T, pow_exponent::Float64=1.0) where {N, T<:NTuple{1,Union{Real,FreeParameter}}} = new(angle, pow_exponent)
end
matrix_rep_raw(g::GPhase{N}) where {N} = Diagonal(SVector{2^N, ComplexF64}(exp(im*g.angle[1]) for _ in 1:2^N))
qubit_count(g::GPhase{N}) where {N} = N
StructTypes.constructfrom(::Type{GPhase}, nt::Quasar.CircuitInstruction) = GPhase{length(nt.targets)}(only(nt.arguments), nt.exponent)

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
        apply_gate!(::Val{$V}, g::GPhase{N}, state_vec::AbstractStateVector{T}, ts::Vararg{Int,N}) where {T<:Complex, N} = apply_gate!($f(exp(im*g.angle[1]*g.pow_exponent)), state_vec, ts...)
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
qubit_count(c::Control{GPhase{N}, B}) where {N, B} = N
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
function matrix_rep_raw(c::Control{GPhase{N}, B}) where {N, B}
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
