"""
    DoubleExcitation(ϕ)

Generate the matrix representation of the [DoubleExcitation](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitation.html) gate.

This gate performs an SO(2) rotation in the subspace {|1100⟩, |0011⟩}, transforming the states as follows:

- |0011⟩ ⟼ cos(ϕ/2)|0011⟩ + sin(ϕ/2)|1100⟩
- |1100⟩ ⟼ cos(ϕ/2)|1100⟩ - sin(ϕ/2)|0011⟩

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = DoubleExcitation(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
julia> eq2 = m * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0];
julia> eq1 ==  [0, 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, 0]
true
julia> eq2 ==  [0, 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, 0]
true
```
"""
struct DoubleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitation(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitation}) = "G2(ang)"
Braket.qubit_count(::Type{DoubleExcitation}) = 4
Base.inv(g::DoubleExcitation) = DoubleExcitation(-g.angle[1])
Base.:^(g::DoubleExcitation, power::Integer) = power == -1 ? inv(g) : (power == 0 ? DoubleExcitation((0.0,)) : DoubleExcitation((g.angle[1] * power,)))
function matrix_rep(g::DoubleExcitation)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)

    mat = diagm(ones(ComplexF64, 16))
    mat[4, 4]   = cosϕ
    mat[13, 13] = cosϕ
    mat[4, 13]  = -sinϕ
    mat[13, 4]  = sinϕ
    return SMatrix{16,16,ComplexF64}(mat)
end

"""
    DoubleExcitationPlus(ϕ)

Generate the matrix representation of the [DoubleExcitationPlus](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitationPlus.html) gate.

This gate performs an SO(2) rotation in the subspace {|1100⟩, |0011⟩} with a phase-shift on other states:

- |0011⟩ ⟼ cos(ϕ/2)|0011⟩ - sin(ϕ/2)|1100⟩
- |1100⟩ ⟼ cos(ϕ/2)|1100⟩ + sin(ϕ/2)|0011⟩
- |x⟩ ⟼ e^{iϕ/2}|x⟩   for all other basis states |x⟩

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = DoubleExcitationPlus(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
julia> eq2 = m * [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];
julia> eq1 ==  [exp(im*ϕ/2), 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, exp(im*ϕ/2)]
true
julia> eq2 ==  [exp(im*ϕ/2), 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, exp(im*ϕ/2)]
true
```

"""
struct DoubleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitationPlus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitationPlus}) = "G2+(ang)"
Braket.qubit_count(::Type{DoubleExcitationPlus}) = 4
Base.inv(g::DoubleExcitationPlus) = DoubleExcitationPlus(-g.angle[1])
Base.:^(g::DoubleExcitationPlus, power::Integer) = power == -1 ? inv(g) : (power == 0 ? DoubleExcitationPlus((0.0,)) : DoubleExcitationPlus((g.angle[1] * power,)))
function matrix_rep(g::DoubleExcitationPlus)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    eiϕ2 = exp(im * g.angle[1] / 2.0)
    mat = diagm(eiϕ2 * ones(ComplexF64, 16))

    mat[4, :] .= 0
    mat[:, 4] .= 0
    mat[13, :] .= 0
    mat[:, 13] .= 0
    # Apply phase-shift to states outside rotation subspace
    mat[4, 4]   = cosϕ
    mat[13, 13] = cosϕ
    mat[4, 13]  = -sinϕ
    mat[13, 4]  = sinϕ
    return SMatrix{16, 16, ComplexF64}(mat)
end
"""
    DoubleExcitationMinus(ϕ)

Generate the matrix representation of the [DoubleExcitationMinus](https://docs.pennylane.ai/en/stable/code/api/pennylane.DoubleExcitationMinus.html) gate.

This gate performs an SO(2) rotation in the subspace {|1100⟩, |0011⟩} with a phase-shift on other states:

- |0011⟩ ⟼ cos(ϕ/2)|0011⟩ - sin(ϕ/2)|1100⟩
- |1100⟩ ⟼ cos(ϕ/2)|1100⟩ + sin(ϕ/2)|0011⟩
- |x⟩ ⟼ e^{-iϕ/2}|x⟩   for all other basis states |x⟩

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = DoubleExcitationMinus(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];
julia> eq2 = m * [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1];
julia> eq1 ==  [exp(-im*ϕ/2), 0, 0, cos(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, sin(ϕ/2), 0, 0, exp(-im*ϕ/2)]
true
julia> eq2 ==  [exp(-im*ϕ/2), 0, 0, -sin(ϕ/2), 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2), 0, 0, exp(-im*ϕ/2)]
true
```

"""
struct DoubleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitationMinus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitationMinus}) = "G2-(ang)"
Braket.qubit_count(::Type{DoubleExcitationMinus}) = 4
Base.inv(g::DoubleExcitationMinus) = DoubleExcitationMinus(-g.angle[1])
Base.:^(g::DoubleExcitationMinus, power::Integer) = power == -1 ? inv(g) : (power == 0 ? DoubleExcitationMinus((0.0,)) : DoubleExcitationMinus((g.angle[1] * power,)))
function matrix_rep(g::DoubleExcitationMinus)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    eiϕ2 = exp(-im * g.angle[1] / 2.0)
    mat = diagm(eiϕ2 * ones(ComplexF64, 16))
    
    mat[4, :] .= 0
    mat[:, 4] .= 0
    mat[13, :] .= 0
    mat[:, 13] .= 0
    # Apply phase-shift to states outside rotation subspace
    mat[4, 4]   = cosϕ
    mat[13, 13] = cosϕ
    mat[4, 13]  = -sinϕ
    mat[13, 4]  = sinϕ
    return SMatrix{16, 16, ComplexF64}(mat)
end

"""
    SingleExcitation(ϕ)

Generate the matrix representation of the [SingleExcitation](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitation.html) gate. 

This gate performs a rotation in the subspace {|01⟩, |10⟩}.

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = SingleExcitation(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [0, 1, 0, 0];
julia> eq2 = m * [0, 0, 1, 0];
julia> eq1 == [0, cos(ϕ/2), - sin(ϕ/2), 0]
true 
julia> eq2 == [0, sin(ϕ/2), cos(ϕ/2), 0]
true
```
"""
struct SingleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitation(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitation}) = "G(ang)"
Braket.qubit_count(::Type{SingleExcitation}) = 2
Base.inv(g::SingleExcitation) = SingleExcitation(-g.angle[1])
Base.:^(g::SingleExcitation, power::Integer) = power == -1 ? inv(g) : (power == 0 ? SingleExcitation((0.0,)) : SingleExcitation((g.angle[1] * power,)))
function matrix_rep(g::SingleExcitation)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([1.0 0 0 0; 0 cosϕ sinϕ 0; 0 -sinϕ cosϕ 0; 0 0 0 1.0])
end

"""
    SingleExcitationPlus(ϕ)

Generate the matrix representation of the [SingleExcitationPlus](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitationPlus.html) gate.

This gate performs a rotation in the subspace {|01⟩, |10⟩} with a phase-shift.

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = SingleExcitationPlus(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [1, 1, 0, 0];
julia> eq2 = m * [1, 0, 1, 0];
julia> eq1 == [exp(im*ϕ/2), cos(ϕ/2), - sin(ϕ/2), 0]
true 
julia> eq2 == [exp(im*ϕ/2), sin(ϕ/2), cos(ϕ/2), 0]
true
```
"""
struct SingleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitationPlus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitationPlus}) = "G+(ang)"
Braket.qubit_count(::Type{SingleExcitationPlus}) = 2
Base.inv(g::SingleExcitationPlus) = SingleExcitationPlus(-g.angle[1])
Base.:^(g::SingleExcitationPlus, power::Integer) = power == -1 ? inv(g) : (power == 0 ? SingleExcitationPlus((0.0,)) : SingleExcitationPlus((g.angle[1] * power,)))
function matrix_rep(g::SingleExcitationPlus)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    eiϕ2 = exp(im * g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([eiϕ2 0 0 0; 0 cosϕ sinϕ 0; 0 -sinϕ cosϕ 0; 0 0 0 eiϕ2])
end

"""
    SingleExcitationMinus(ϕ)

Generate the matrix representation of the [SingleExcitationMinus](https://docs.pennylane.ai/en/stable/code/api/pennylane.SingleExcitationMinus.html) gate.

This gate performs a rotation in the subspace {|01⟩, |10⟩} with a phase-shift.

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix = SingleExcitationMinus(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [1, 1, 0, 0];
julia> eq2 = m * [1, 0, 1, 0];
julia> eq1 == [exp(-im*ϕ/2), cos(ϕ/2), - sin(ϕ/2), 0]
true 
julia> eq2 == [exp(-im*ϕ/2), sin(ϕ/2), cos(ϕ/2), 0]
true
```

"""
struct SingleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitationMinus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitationMinus}) = "G-(ang)"
Braket.qubit_count(::Type{SingleExcitationMinus}) = 2
Base.inv(g::SingleExcitationMinus) = SingleExcitationMinus(-g.angle[1])
Base.:^(g::SingleExcitationMinus, power::Integer) = power == -1 ? inv(g) : (power == 0 ? SingleExcitationMinus((0.0,)) : SingleExcitationMinus((g.angle[1] * power,)))
function matrix_rep(g::SingleExcitationMinus)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    eiϕ2 = exp(-im * g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([eiϕ2 0 0 0; 0 cosϕ sinϕ 0; 0 -sinϕ cosϕ 0; 0 0 0 eiϕ2])
end

"""
    FermionicSWAP(ϕ)

Generate the matrix representation of the [FermionicSWAP](https://docs.pennylane.ai/en/stable/code/api/pennylane.FermionicSWAP.html) gate.

This gate performs a rotation in adjacent fermionic modes under the Jordan-Wigner mapping, transforming states as follows:

- |00⟩ ⟼ |00⟩
- |01⟩ ⟼ e^{iϕ/2}cos(ϕ/2)|01⟩ - ie^{iϕ/2}sin(ϕ/2)|10⟩
- |10⟩ ⟼ -ie^{iϕ/2}sin(ϕ/2)|01⟩ + e^{iϕ/2}cos(ϕ/2)|10⟩
- |11⟩ ⟼ e^{iϕ}|11⟩

# Examples

```jldoctest
julia> ϕ = 3.56;
julia> gate_matrix =  FermionicSWAP(ϕ);
julia> m  = matrix_rep(gate_matrix);
julia> eq1 = m * [0, 0, 0, 0];
julia> eq2 = m * [0, 1, 0, 0];
julia> eq3 = m * [0, 0, 1, 0];
julia> eq4 = m * [0, 0, 0, 1];
julia> eq1 == [0, 0, 0, 0]
true 
julia> eq2 == [0, exp(im*ϕ/2.0)*cos(ϕ / 2.0), - im*exp(im*ϕ/2.0)*sin(ϕ/2.0), 0]
true
julia> eq3 == [0, - im*exp(im*ϕ/2.0)*sin(ϕ/2.0), exp(im*ϕ/2.0)*cos(ϕ/2.0), 0]
true
julia> eq4 == [0, 0, 0, exp(im * ϕ)]
true
```
"""
struct FermionicSWAP <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    FermionicSWAP(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{FermionicSWAP}) = "FSWAP(ang)"
Braket.qubit_count(::Type{FermionicSWAP}) = 2
Base.inv(g::FermionicSWAP) = FermionicSWAP(-g.angle[1])
Base.:^(g::FermionicSWAP, power::Integer) = power == -1 ? inv(g) : (power == 0 ? FermionicSWAP((0.0,)) : FermionicSWAP((g.angle[1] * power,)))
function matrix_rep(g::FermionicSWAP)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    eiϕ2 = exp(im * g.angle[1] / 2.0)
    eiϕ = exp(im * g.angle[1])
    ieiϕ2 = im * eiϕ2

    return SMatrix{4,4,ComplexF64}([1.0 0 0 0; 0 eiϕ2 * cosϕ -ieiϕ2 * sinϕ 0; 0 -ieiϕ2 * sinϕ eiϕ2 * cosϕ 0; 0 0 0 eiϕ])
end

struct MultiRZ <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    MultiRZ(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} = new(angle)
end
Braket.chars(::Type{MultiRZ}) = "MultiRZ(ang)"
Braket.qubit_count(g::MultiRZ) = 1
Base.inv(g::MultiRZ) = MultiRZ(-g.angle[1])

struct U <: AngledGate{3}
    angle::NTuple{3,Union{Float64,FreeParameter}}
    U(angle::T) where {T<:NTuple{3,Union{Float64,FreeParameter}}} =
        new(angle)
end
U(θ::T, ϕ::T, λ::T) where {T<:Union{Float64,FreeParameter}} = U((θ, ϕ, λ))
function matrix_rep(g::U)
    θ, ϕ, λ = g.angle
    return SMatrix{2, 2, ComplexF64}([cos(θ/2) -exp(im*λ)*sin(θ/2); exp(im*ϕ)*sin(θ/2) exp(im*(λ+ϕ))*cos(θ/2)])
end
Braket.chars(::Type{U}) = "U(ang)"
Braket.qubit_count(g::U) = 1
Base.inv(g::U) = U(-g.angle[1], -g.angle[3], -g.angle[2])

struct MultiQubitPhaseShift{N} <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    MultiQubitPhaseShift{N}(angle::NTuple{1,<:Union{Float64,FreeParameter}}) where {N} = new(angle)
end
matrix_rep(g::MultiQubitPhaseShift{N}) where {N} = Diagonal(SVector{2^N, ComplexF64}(exp(im*g.angle[1]) for _ in 1:2^N))
Braket.chars(::Type{MultiQubitPhaseShift}) = "GPhase(ang)"
Braket.qubit_count(g::MultiQubitPhaseShift{N}) where {N} = N
Base.inv(g::MultiQubitPhaseShift{N}) where {N} = MultiQubitPhaseShift{N}((-g.angle[1],))

function apply_gate!(
    factor::ComplexF64,
    diagonal::Braket.PauliEigenvalues{N},
    state_vec::StateVector{T},
    ts::Vararg{Int,N},
) where {T<:Complex,N}
    g_mat  = Diagonal(SVector{2^N,ComplexF64}(exp(factor * diagonal[i]) for i = 1:2^N))
    apply_gate!(g_mat, state_vec, ts...)
end
function apply_gate!(
    factor::ComplexF64,
    state_vec::StateVector{T},
    ts::Vararg{Int,N},
) where {T<:Complex,N}
    g_mat = Diagonal(SVector{2^N,ComplexF64}(factor for i = 1:2^N))
    apply_gate!(g_mat, state_vec, ts...)
end
for (V, f) in ((true, :conj), (false, :identity))
    @eval begin
        apply_gate!(
            ::Val{$V},
            g::MultiRZ,
            state_vec::StateVector{T},
            t::Int,
        ) where {T<:Complex} = apply_gate!(Val($V), Rz(g.angle), state_vec, t)
        apply_gate!(
            ::Val{$V},
            g::MultiRZ,
            state_vec::StateVector{T},
            t1::Int,
            t2::Int,
        ) where {T<:Complex} = apply_gate!(Val($V), ZZ(g.angle), state_vec, t1, t2)
        apply_gate!(
            ::Val{$V},
            g::MultiRZ,
            state_vec::StateVector{T},
            ts::Vararg{Int,N},
        ) where {T<:Complex,N} = apply_gate!($f(-im * g.angle[1] / 2.0), Braket.PauliEigenvalues(Val(N)), state_vec, ts...)
        apply_gate!(::Val{$V}, g::MultiQubitPhaseShift{N}, state_vec::StateVector{T}, ts::Vararg{Int,N}) where {T<:Complex, N} = apply_gate!($f(exp(im*g.angle[1])), state_vec, ts...)
    end
end

function apply_gate!(
    g::DoubleExcitation,
    state_vec::StateVector{T},
    t1::Int,
    t2::Int,
    t3::Int,
    t4::Int,
) where {T<:Complex}
    n_amps, endian_ts = get_amps_and_qubits(state_vec, t1, t2, t3, t4)
    ordered_ts = sort(collect(endian_ts))
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    e_t1, e_t2, e_t3, e_t4 = endian_ts
    Threads.@threads for ix = 0:div(n_amps, 2^4)-1
        padded_ix = pad_bits(ix, ordered_ts)
        i0011 = flip_bits(padded_ix, (e_t3, e_t4)) + 1
        i1100 = flip_bits(padded_ix, (e_t1, e_t2)) + 1
        @inbounds begin
            amp0011 = state_vec[i0011]
            amp1100 = state_vec[i1100]
            state_vec[i0011] = cosϕ * amp0011 - sinϕ * amp1100
            state_vec[i1100] = sinϕ * amp0011 + cosϕ * amp1100
        end
    end
    return
end
apply_gate!(::Val{true}, g::DoubleExcitation, state_vec::StateVector{T}, t1::Int, t2::Int, t3::Int, t4::Int) where {T<:Complex} = apply_gate!(g, state_vec, t1, t2, t3, t4)
apply_gate!(::Val{false}, g::DoubleExcitation, state_vec::StateVector{T}, t1::Int, t2::Int, t3::Int, t4::Int) where {T<:Complex} = apply_gate!(g, state_vec, t1, t2, t3, t4)

struct Control{G<:Gate, B} <: Gate
    g::G
    bitvals::NTuple{B, Int}
    Control{G, B}(g::G, bitvals::NTuple{B, Int}) where {B, G} = new(g, bitvals)
end
Control(g::G, bitvals::NTuple{B, Int}) where {G<:Gate, B} = Control{G, B}(g, bitvals)
Control(g::Control{G, BC}, bitvals::NTuple{B, Int}) where {G<:Gate, BC, B} = Control(g.g, (g.bitvals..., bitvals...))
Control(g::G, bitvals::NTuple{0, Int}) where {G<:Gate} = g
Braket.qubit_count(c::Control{G, B}) where {G<:Gate, B} = qubit_count(c.g) + B
Braket.qubit_count(c::Control{MultiQubitPhaseShift{N}, B}) where {N, B} = N
function matrix_rep(c::Control{G, B}) where {G<:Gate, B}
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
function matrix_rep(c::Control{MultiQubitPhaseShift{N}, B}) where {N, B}
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

const custom_gates = (doubleexcitation=DoubleExcitation, singleexcitation=SingleExcitation, multirz=MultiRZ, u=U,) 

