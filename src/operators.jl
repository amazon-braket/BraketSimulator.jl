"""
    Operator

Abstract type representing operations that can be applied to a [`Circuit`](@ref).
Subtypes include [`Gate`](@ref), [`Noise`](@ref), [`Observable`](@ref BraketSimulator.Observables.Observable).
"""
abstract type Operator end

abstract type Parametrizable end
struct Parametrized end
struct NonParametrized end

"""
    QuantumOperator < Operator

Abstract type representing *quantum* operations that can be applied to a [`Circuit`](@ref).
Subtypes include [`Gate`](@ref) and [`Noise`](@ref).
"""
abstract type QuantumOperator <: Operator end
StructTypes.StructType(::Type{QuantumOperator}) = StructTypes.AbstractType()
StructTypes.subtypes(::Type{QuantumOperator}) = (h=H, i=I, x=X, y=Y, z=Z, s=S, si=Si, t=T, ti=Ti, v=V, vi=Vi, cnot=CNot, swap=Swap, iswap=ISwap, cv=CV, cy=CY, cz=CZ, ecr=ECR, ccnot=CCNot, cswap=CSwap, unitary=Unitary, rx=Rx, ry=Ry, rz=Rz, phaseshift=PhaseShift, pswap=PSwap, xy=XY, cphaseshift=CPhaseShift, cphaseshift00=CPhaseShift00, cphaseshift01=CPhaseShift01, cphaseshift10=CPhaseShift10, xx=XX, yy=YY, zz=ZZ, gpi=GPi, gpi2=GPi2, ms=MS, prx=PRx, u=U, gphase=GPhase, kraus=Kraus, bit_flip=BitFlip, phase_flip=PhaseFlip, pauli_channel=PauliChannel, amplitude_damping=AmplitudeDamping, phase_damping=PhaseDamping, depolarizing=Depolarizing, two_qubit_dephasing=TwoQubitDephasing, two_qubit_depolarizing=TwoQubitDepolarizing, generalized_amplitude_damping=GeneralizedAmplitudeDamping, multi_qubit_pauli_channel=MultiQubitPauliChannel, measure=Measure, reset=Reset, barrier=Barrier, delay=Delay)
parameters(::QuantumOperator) = FreeParameter[]

qubit_count(o::Matrix) = Int(log2(size(o, 1))) 

struct PauliEigenvalues{N}
    coeff::Float64
    PauliEigenvalues{N}(coeff::Float64=1.0) where {N} = new(coeff)
end
PauliEigenvalues(::Val{N}, coeff::Float64=1.0) where {N} = PauliEigenvalues{N}(coeff)
Base.length(::PauliEigenvalues{N}) where {N} = 2^N
Base.iterate(p::PauliEigenvalues{N}, ix::Int=1) where {N} = ix <= length(p) ? (p[ix], ix+1) : nothing
Base.getindex(p::PauliEigenvalues{1}, i::Int)::Float64 = getindex((p.coeff, -p.coeff), i)
function Base.getindex(p::PauliEigenvalues{N}, i::Int)::Float64 where N
    i_block = div(i-1, 2)
    split = div(2^(N-1)-1, 2)
    if N < 5
        is_front = !isodd(mod(i-1, 2))
        ev = is_front ? p.coeff : -p.coeff
        mi = mod(i_block, 2)
        di = div(i_block, 2)
        if i_block <= split
            return isodd(mi) ⊻ isodd(di) ? -ev : ev
        else
            mi = mod(i_block - split - 1, 2)
            di = div(i_block - split - 1, 2)
            return isodd(mi) ⊻ isodd(di) ? ev : -ev
        end
    else
        new_i = i > 2^(N-1) ? i - 2^(N-1) : i
	    one_down_pe = PauliEigenvalues(Val(N-1))
	    one_down = one_down_pe[new_i]
        return i_block <= split ? one_down : -one_down
    end
end
Base.getindex(p::PauliEigenvalues{N}, ix::Vector{Int}) where {N} = [p[i] for i in ix]

"""
    Measure(index, result) <: QuantumOperator

Represents a measurement operation on targeted qubit, stored in the classical register at `index`.
The `result` field represents the measurement outcome (0 or 1) and allows the operator to act as a projector
into the measured subspace.
"""
struct Measure <: QuantumOperator
    index::Int
    result::Int
end
Measure() = Measure(-1, -1)
Measure(index::Int) = Measure(index, -1)
StructTypes.constructfrom(::Type{Measure}, nt) = Measure()

"""
    Reset() <: QuantumOperator

Represents an active reset operation on targeted qubit.
For now, this is a no-op.
"""
struct Reset <: QuantumOperator end
StructTypes.constructfrom(::Type{Reset}, nt) = Reset()

"""
    Barrier() <: QuantumOperator

Represents a barrier operation on targeted qubit.
For now, this is a no-op.
"""
struct Barrier <: QuantumOperator end
StructTypes.constructfrom(::Type{Barrier}, nt) = Barrier()

"""
    Delay(duration::Time) <: QuantumOperator

Represents a delay operation for `duration` on targeted qubit.
For now, this is a no-op.
"""
struct Delay <: QuantumOperator
    duration::Dates.Period
end
StructTypes.constructfrom(::Type{Delay}, nt) = Delay(only(nt.arguments))

for T in (:Barrier, :Reset, :Delay, :Measure)
    @eval begin
        qubit_count(::Type{$T}) = 1
        qubit_count(o::$T) = qubit_count($T)
        Parametrizable(::$T) = NonParametrized()
    end
end
