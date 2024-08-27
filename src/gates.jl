"""
    Gate <: QuantumOperator

Abstract type representing a quantum gate.
"""
abstract type Gate <: QuantumOperator end

"""
    AngledGate{NA} <: Gate

Parametric type representing a quantum gate with `NA` `angle` parameters.
"""
abstract type AngledGate{NA} <: Gate end
n_angles(::Type{<:Gate}) = 0
n_angles(::Type{<:AngledGate{N}}) where {N} = N
n_angles(g::G) where {G<:Gate} = n_angles(G)

for gate_def in (
                 (:Rx, :1, :1, :0, :1),
                 (:Ry, :1, :1, :0, :1),
                 (:Rz, :1, :1, :0, :1),
                 (:PSwap, :1, :2, :0, :2),
                 (:PhaseShift, :1, :1, :0, :1),
                 (:CPhaseShift, :1, :2, :1, :1),
                 (:CPhaseShift00, :1, :2, :1, :1),
                 (:CPhaseShift01, :1, :2, :1, :1),
                 (:CPhaseShift10, :1, :2, :1, :1),
                 (:XX, :1, :2, :0, :2),
                 (:XY, :1, :2, :0, :2),
                 (:YY, :1, :2, :0, :2),
                 (:ZZ, :1, :2, :0, :2),
                 (:GPi, :1, :1, :0, :1),
                 (:GPi2, :1, :1, :0, :1),
                 (:MS, :3, :2, :0, :2),
                 (:U, :3, :1, :0, :1),
                 (:PRx, :2, :1, :0, :1),
                )
    G, n_angle, qc, n_controls, n_targets = gate_def
    @eval begin
        @doc """
            $($G) <: AngledGate{$($n_angle)}
            $($G)(angles) -> $($G)
        
        $($G) gate.
        """
        mutable struct $G <: AngledGate{$n_angle}
            angle::NTuple{$n_angle, Union{Real, FreeParameter}}
            pow_exponent::Float64
            $G(angle::T, pow_exponent::Float64=1.0) where {T<:NTuple{$n_angle, Union{Real, FreeParameter}}} = new(angle, pow_exponent)
        end
        qubit_count(::Type{$G}) = $qc
    end
end
Base.:(==)(g1::G, g2::G) where {G<:AngledGate} = g1.pow_exponent == g2.pow_exponent && g1.angle == g2.angle

for gate_def in (
                 (:H, :1, :0, :1),
                 (:I, :1, :0, :1),
                 (:X, :1, :0, :1),
                 (:Y, :1, :0, :1),
                 (:Z, :1, :0, :1),
                 (:S, :1, :0, :1),
                 (:Si, :1, :0, :1),
                 (:T, :1, :0, :1),
                 (:Ti, :1, :0, :1),
                 (:V, :1, :0, :1),
                 (:Vi, :1, :0, :1),
                 (:CNot, :2, :1, :1),
                 (:Swap, :2, :0, :2),
                 (:ISwap, :2, :0, :2),
                 (:CV, :2, :1, :1),
                 (:CY, :2, :1, :1),
                 (:CZ, :2, :1, :1),
                 (:ECR, :2, :0, :2),
                 (:CCNot, :3, :2, :1),
                 (:CSwap, :3, :1, :2),
                )
    G, qc, n_controls, n_targets = gate_def
    @eval begin
        @doc """
            $($G) <: Gate
            $($G)() -> $($G)
        
        $($G) gate.
        """
        mutable struct $G <: Gate
            pow_exponent::Float64
            $G(pow_exponent::Float64=1.0) = new(pow_exponent)
        end
        qubit_count(::Type{$G}) = $qc
    end
end
(::Type{G})(angle::T, pow_exponent=1.0) where {G<:AngledGate{1}, T<:Union{Real, FreeParameter}} = G((angle,), Float64(pow_exponent))
(::Type{G})(angle1::T1, angle2::T2, pow_exponent=1.0) where {T1<:Union{Real, FreeParameter}, T2<:Union{Real, FreeParameter}, G<:AngledGate{2}} = G((angle1, angle2,), Float64(pow_exponent))
(::Type{G})(angle1::T1, angle2::T2, angle3::T3, pow_exponent=1.0) where {T1<:Union{Real, FreeParameter}, T2<:Union{Real, FreeParameter}, T3<:Union{Real, FreeParameter}, G<:AngledGate{3}} = G((angle1, angle2, angle3,), Float64(pow_exponent))
Base.:(==)(g1::G, g2::G) where {G<:Gate} = g1.pow_exponent == g2.pow_exponent
qubit_count(g::G) where {G<:Gate}  = qubit_count(G)
angles(g::G) where {G<:Gate}       = ()
angles(g::AngledGate{N}) where {N} = g.angle

"""
    Unitary <: Gate
    Unitary(matrix::Matrix{ComplexF64}) -> Unitary

Arbitrary unitary gate.
"""
mutable struct Unitary <: Gate
    matrix::Matrix{ComplexF64}
    pow_exponent::Float64
    Unitary(matrix::Matrix{<:Number}, pow_exponent=1.0) = new(ComplexF64.(matrix), Float64(pow_exponent))
end
Base.:(==)(u1::Unitary, u2::Unitary) = u1.matrix == u2.matrix && u1.pow_exponent == u2.pow_exponent
qubit_count(g::Unitary) = convert(Int, log2(size(g.matrix, 1)))

Parametrizable(g::AngledGate) = Parametrized()
Parametrizable(g::Gate)       = NonParametrized()
parameters(g::AngledGate)     = collect(filter(a->a isa FreeParameter, angles(g)))
# nosemgrep
function bind_value!(::Parametrized, g::G, params::Dict{Symbol, <:Real}) where {G<:AngledGate}
    new_angles = map(angles(g)) do angle
        angle isa FreeParameter || return angle
        return get(params, angle.name, angle)
    end
    return G(new_angles...)
end
Base.copy(g::G) where {G<:Gate} = G((copy(getproperty(g, fn)) for fn in fieldnames(G))...)
Base.copy(g::G) where {G<:AngledGate} = G(angles(g), g.pow_exponent)
