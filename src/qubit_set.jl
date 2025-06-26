using OrderedCollections

mutable struct Qubit <: Integer
    index::Int
    name::String
    measured::Bool
    Qubit(q::Integer) = new(q, "q" * string(q), false)
    Qubit(q::AbstractFloat) = new(Int(q), "q" * string(Int(q)), false)
    Qubit(q::BigFloat) = new(Int(q), "q" * string(Int(q)), false)
    Qubit(q::Integer, name::String) = new(q, name, false)
    Qubit(q::Integer, name::String, measured::Bool) = new(q, name, measured)
end
Qubit(q::Qubit) = new(q.index, q.name, q.measured)
Base.:(==)(q::Qubit, i::T) where {T<:Integer} = q.index==i
Base.:(==)(i::T, q::Qubit) where {T<:Integer} = q.index==i
Base.:(==)(i::BigInt, q::Qubit) = big(q.index)==i
Base.:(==)(q::Qubit, i::BigInt) = big(q.index)==i
Base.:(==)(q1::Qubit, q2::Qubit) = q1.index==q2.index

Base.hash(q::Qubit, h::UInt) = hash(q.index, h)
Base.show(io::IO, q::Qubit) = print(io, "Qubit($(q.name))")
const IntOrQubit    = Union{Int, Qubit}

"""
    QubitSet

An `OrderedSet`-like object which represents the qubits a
[`Circuit`](@ref), [`Instruction`](@ref), or [`Result`](@ref)
acts on and their ordering.

# Examples
```jldoctest
julia> QubitSet([2, 1])
QubitSet with 2 elements:
  2
  1

julia> QubitSet()
QubitSet()

julia> QubitSet(QubitSet(5, 1))
QubitSet with 2 elements:
  5
  1
```
"""
struct QubitSet <: AbstractSet{Int}
    dict::OrderedDict{IntOrQubit, Nothing}
    QubitSet()   = new(OrderedDict{IntOrQubit, Nothing}())
    QubitSet(xs) = (v = foldl(vcat, xs, init=Int[]); return union!(new(OrderedDict{IntOrQubit,Nothing}()), v))
    QubitSet(qs::Vararg{IntOrQubit}) = QubitSet(collect(qs))
    QubitSet(qs::QubitSet) = qs
    QubitSet(::Nothing) = QubitSet()
end
Base.convert(::Type{Vector{Int}}, qs::QubitSet) = convert.(Int, collect(qs))
Base.length(qs::QubitSet)   = length(qs.dict)
Base.lastindex(qs::QubitSet) = length(qs)
Base.isempty(qs::QubitSet)  = isempty(qs.dict)
Base.in(q, qs::QubitSet)    = haskey(qs.dict, q)
Base.push!(qs::QubitSet, q) = (qs.dict[q] = nothing; qs)
Base.copy(qs::QubitSet)     = QubitSet(qs[ii] for ii in 1:length(qs))
Base.popfirst!(qs::QubitSet) = (q = popfirst!(qs.dict); return q[1])
function Base.iterate(qs::QubitSet)::Union{Nothing, Tuple{IntOrQubit, Int}}
    qs.dict.ndel > 0 && OrderedCollections.rehash!(qs.dict)
    length(qs.dict.keys) < 1 && return nothing
    return (qs.dict.keys[1], 2)
end
function Base.iterate(qs::QubitSet, i)::Union{Nothing, Tuple{IntOrQubit, Int}}
    length(qs.dict.keys) < i && return nothing
    return (qs.dict.keys[i], i+1)
end
Base.:(==)(q1::QubitSet, q2::QubitSet) = (length(q1) == length(q2)) && all(q1[ii] == q2[ii] for ii in 1:length(q1))

function Base.getindex(qs::QubitSet, i::Int)
    qs.dict.ndel > 0 && OrderedCollections.rehash!(qs.dict)
    return qs.dict.keys[i]
end
Base.getindex(qs::QubitSet, ui::UnitRange) = QubitSet([qs[ii] for ii in ui])

function Base.intersect(qs1::QubitSet, qs2::QubitSet)
    qs = QubitSet()
    for q in qs1
        (q in qs2) && union!(qs, q)
    end
    return qs
end

function Base.show(io::IO, qs::QubitSet)
    print(io, "QubitSet(")
    q_strs = map(q->sprint(show, q), collect(qs)) 
    print(io, join(q_strs, ", "))
    print(io, ")")
end
Base.convert(::Type{QubitSet}, v::Vector{<:Integer}) = QubitSet(v)
Base.sort(qs::QubitSet; kwargs...) = QubitSet(sort(collect(qs); kwargs...))

const VecOrQubitSet = Union{Vector, QubitSet}
