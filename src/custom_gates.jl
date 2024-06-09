struct DoubleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitation(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitation}) = "G2(ang)"
Braket.qubit_count(::Type{DoubleExcitation}) = 4
Base.inv(g::DoubleExcitation) = DoubleExcitation(-g.angle[1])
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

struct DoubleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitationPlus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitationPlus}) = "G2+(ang)"
Braket.qubit_count(::Type{DoubleExcitationPlus}) = 4
Base.inv(g::DoubleExcitationPlus) = DoubleExcitationPlus(-g.angle[1])
function matrix_rep(g::DoubleExcitationPlus)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    eᵢᵩ₂ = exp(im * g.angle[1] / 2.0)
 
    mat = diagm(fill(eᵢᵩ₂, 16))
    mat[4, 4]   = cosᵩ
    mat[13, 13] = cosᵩ
    mat[4, 13]  = -sinᵩ
    mat[13, 4]  = sinᵩ
    return SMatrix{16,16,ComplexF64}(mat)
end

struct DoubleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    DoubleExcitationMinus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{DoubleExcitationMinus}) = "G2-(ang)"
Braket.qubit_count(::Type{DoubleExcitationMinus}) = 4
Base.inv(g::DoubleExcitationMinus) = DoubleExcitationMinus(-g.angle[1])
function matrix_rep(g::DoubleExcitationMinus)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    eᵢᵩ₂ = exp(-im * g.angle[1] / 2.0)

    mat = diagm(fill(eᵢᵩ₂, 16))
    mat[4, 4]   = cosᵩ
    mat[13, 13] = cosᵩ
    mat[4, 13]  = -sinᵩ
    mat[13, 4]  = sinᵩ
    return SMatrix{16,16,ComplexF64}(mat)
end


struct SingleExcitation <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitation(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitation}) = "G(ang)"
Braket.qubit_count(::Type{SingleExcitation}) = 2
Base.inv(g::SingleExcitation) = SingleExcitation(-g.angle[1])
function matrix_rep(g::SingleExcitation)
    cosϕ = cos(g.angle[1] / 2.0)
    sinϕ = sin(g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([1.0 0 0 0; 0 cosϕ sinϕ 0; 0 -sinϕ cosϕ 0; 0 0 0 1.0])
end

struct SingleExcitationPlus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitationPlus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitationPlus}) = "G+(ang)"
Braket.qubit_count(::Type{SingleExcitationPlus}) = 2
Base.inv(g::SingleExcitationPlus) = SingleExcitationPlus(-g.angle[1])
function matrix_rep(g::SingleExcitationPlus)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    eᵢᵩ₂ = exp(im * g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([eᵢᵩ₂ 0 0 0; 0 cosᵩ sinᵩ 0; 0 -sinᵩ cosᵩ 0; 0 0 0 eᵢᵩ₂])
end

struct SingleExcitationMinus <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    SingleExcitationMinus(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{SingleExcitationMinus}) = "G-(ang)"
Braket.qubit_count(::Type{SingleExcitationMinus}) = 2
Base.inv(g::SingleExcitationMinus) = SingleExcitationMinus(-g.angle[1])
function matrix_rep(g::SingleExcitationMinus)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    eᵢᵩ₂ = exp(-im * g.angle[1] / 2.0)
    return SMatrix{4,4,ComplexF64}([eᵢᵩ₂ 0 0 0; 0 cosᵩ sinᵩ 0; 0 -sinᵩ cosᵩ 0; 0 0 0 eᵢᵩ₂])
end

struct OrbitalRotation <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    OrbitalRotation(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{OrbitalRotation}) = "O(ang)"
Braket.qubit_count(::Type{OrbitalRotation}) = 4
Base.inv(g::OrbitalRotation) = OrbitalRotation(-g.angle[1])
function matrix_rep(g::OrbitalRotation)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    # This transformation applies to the two neighboring spatial orbitals.
    mat = SMatrix{16,16,ComplexF64}(I) # Identity matrix
    # Define the subspace transformation
    mat[1, 1] = cosᵩ
    mat[2, 2] = cosᵩ
    mat[1, 2] = -sinᵩ
    mat[2, 1] = sinᵩ
    mat[15, 15] = cosᵩ
    mat[16, 16] = cosᵩ
    mat[15, 16] = -sinᵩ
    mat[16, 15] = sinᵩ
    return mat
end

struct FermionicSWAP <: AngledGate{1}
    angle::NTuple{1,Union{Float64,FreeParameter}}
    FermionicSWAP(angle::T) where {T<:NTuple{1,Union{Float64,FreeParameter}}} =
        new(angle)
end
Braket.chars(::Type{FermionicSWAP}) = "FSWAP(ang)"
Braket.qubit_count(::Type{FermionicSWAP}) = 2
Base.inv(g::FermionicSWAP) = FermionicSWAP(-g.angle[1])
function matrix_rep(g::FermionicSWAP)
    cosᵩ = cos(g.angle[1] / 2.0)
    sinᵩ = sin(g.angle[1] / 2.0)
    eᵢᵩ₂ = exp(im * g.angle[1] / 2.0)
    eᵢᵩ = exp(im * g.angle[1])
    ieiϕ2 = im * eiϕ2
   
    return SMatrix{4,4,ComplexF64}([
        1 0 0 0 0;
        0 eᵢᵩ₂*cosᵩ -im*eᵢᵩ₂ * sinᵩ 0;
        0 -im*eᵢᵩ₂ * sinᵩ eᵢᵩ₂ * cosᵩ 0;
        0 0 0 eᵢᵩ
    ])
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

