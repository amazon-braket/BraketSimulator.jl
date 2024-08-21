"""
    Operator

Abstract type representing operations that can be applied to a [`Circuit`](@ref).
Subtypes include [`Gate`](@ref), [`Noise`](@ref), [`Observable`](@ref BraketSimulator.Observables.Observable).
"""
abstract type Operator end

"""
    QuantumOperator < Operator

Abstract type representing *quantum* operations that can be applied to a [`Circuit`](@ref).
Subtypes include [`Gate`](@ref) and [`Noise`](@ref).
"""
abstract type QuantumOperator <: Operator end

abstract type Parametrizable end
struct Parametrized end 
struct NonParametrized end 

struct PauliEigenvalues{N}
    coeff::Float64
    PauliEigenvalues{N}(coeff::Float64=1.0) where {N} = new(coeff)
end
PauliEigenvalues(::Val{N}, coeff::Float64=1.0) where {N} = PauliEigenvalues{N}(coeff)
Base.length(p::PauliEigenvalues{N}) where {N} = 2^N
function Base.iterate(p::PauliEigenvalues{N}, ix::Int=1) where {N}
    return ix <= length(p) ? (p[ix], ix+1) : nothing
end

Base.getindex(p::PauliEigenvalues{1}, i::Int)::Float64 = getindex((p.coeff, -p.coeff), i)
function Base.getindex(p::PauliEigenvalues{N}, i::Int)::Float64 where N
    i_block = div(i-1, 2)
    split = div(2^(N-1)-1, 2)
    if N < 5
        total_evs = 2^N
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
    Measure(index) <: QuantumOperator

Represents a measurement operation on targeted qubit, stored in the classical register at `index`.
"""
struct Measure <: QuantumOperator
    index::Int
end
Measure() = Measure(-1)
Parametrizable(m::Measure) = NonParametrized()
qubit_count(::Type{Measure}) = 1
qubit_count(m::Measure) = qubit_count(Measure)

"""
    Reset(index) <: QuantumOperator

Represents an active reset operation on targeted qubit, stored in the classical register at `index`.
For now, this is a no-op.
"""
struct Reset <: QuantumOperator
    index::Int
end
Reset() = Reset(-1)
Parametrizable(m::Reset) = NonParametrized()
qubit_count(::Type{Reset}) = 1
qubit_count(m::Reset) = qubit_count(Reset)

"""
    Barrier(index) <: QuantumOperator

Represents a barrier operation on targeted qubit, stored in the classical register at `index`.
For now, this is a no-op.
"""
struct Barrier <: QuantumOperator
    index::Int
end
Barrier() = Barrier(-1)
Parametrizable(m::Barrier) = NonParametrized()
qubit_count(::Type{Barrier}) = 1
qubit_count(m::Barrier) = qubit_count(Barrier)

"""
    Delay(index, duration::Time) <: QuantumOperator

Represents a delay operation for `duration` on targeted qubit,
stored in the classical register at `index`.
For now, this is a no-op.
"""
struct Delay <: QuantumOperator
    index::Int
    duration::Dates.Period
end
Delay(duration::Dates.Period) = Delay(-1, duration)
Parametrizable(m::Delay) = NonParametrized()
qubit_count(::Type{Delay}) = 1
qubit_count(m::Delay) = qubit_count(Delay)
