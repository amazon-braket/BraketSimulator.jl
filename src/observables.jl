module Observables

using StructTypes, LinearAlgebra

import ..BraketSimulator: qubit_count, complex_matrix_to_ir, complex_matrix_from_ir, Operator, QubitSet, Qubit, IntOrQubit, IRObservable, PauliEigenvalues
export Observable, TensorProduct, HermitianObservable, Sum

"""
    Observable <: Operator

Abstract type representing an observable to be measured. All `Observable`s
have `eigvals` defined.

See also: [`H`](@ref), [`I`](@ref), [`X`](@ref), [`Y`](@ref), [`Z`](@ref), [`TensorProduct`](@ref), [`HermitianObservable`](@ref).
"""
abstract type Observable <: Operator end
abstract type NonCompositeObservable <: Observable end
abstract type StandardObservable <: NonCompositeObservable end

LinearAlgebra.ishermitian(o::Observable) = true

coef(o::O) where {O<:Observable} = o.coefficient

for (typ, label, super) in ((:H, "h", :StandardObservable), (:X, "x", :StandardObservable), (:Y, "y", :StandardObservable), (:Z, "z", :StandardObservable), (:I, "i", :NonCompositeObservable))
    @eval begin
        @doc """
            $($typ) <: Observable
            $($typ)([coeff::Float64]) -> $($typ)

        Struct representing a `$($typ)` observable in a measurement. The observable
        may be scaled by `coeff`.
        """
        struct $typ <: $super
            coefficient::Float64
            $typ(coef::Float64=1.0) = new(coef)
        end
        StructTypes.lower(x::$typ) =  Union{String, Vector{Vector{Vector{Float64}}}}[$label]
        qubit_count(::Type{$typ}) = 1
        qubit_count(o::$typ) = qubit_count($typ)
        unscaled(o::$typ) = $typ()
        Base.copy(o::$typ) = $typ(o.coefficient)
        Base.:(*)(o::$typ, n::Real) = $typ(Float64(n*o.coefficient))
        Base.:(==)(o1::$typ, o2::$typ) = (o1.coefficient ≈ o2.coefficient)
        Base.show(io::IO, o::$typ) = print(io, (isone(o.coefficient) ? "" : string(o.coefficient) * " * ") * uppercase($label))
    end
end
LinearAlgebra.eigvals(o::I) = [o.coefficient, o.coefficient]
LinearAlgebra.eigvals(o::StandardObservable) = [o.coefficient, -o.coefficient]

"""
    HermitianObservable <: Observable
    HermitianObservable(matrix::Matrix) -> HermitianObservable

Struct representing an observable of an arbitrary complex Hermitian matrix.

# Examples
```jldoctest
julia> ho = BraketSimulator.Observables.HermitianObservable([0 1; 1 0])
HermitianObservable((2, 2))
```
"""
struct HermitianObservable <: NonCompositeObservable
    matrix::Matrix{<:Complex}
    coefficient::Float64
    function HermitianObservable(mat::Matrix{<:Number})
        ishermitian(mat) || throw(ArgumentError("input matrix to HermitianObservable must be Hermitian."))
        new(complex(mat), 1.0)
    end
end
HermitianObservable(v::Vector{Vector{Vector{T}}}) where {T<:Number} = HermitianObservable(complex_matrix_from_ir(v))
Base.copy(o::HermitianObservable) = HermitianObservable(copy(o.matrix))
StructTypes.lower(x::HermitianObservable) = Union{String, Vector{Vector{Vector{Float64}}}}[complex_matrix_to_ir(ComplexF64.(x.matrix))]
Base.:(==)(h1::HermitianObservable, h2::HermitianObservable) = (size(h1.matrix) == size(h2.matrix) && h1.matrix ≈ h2.matrix)
qubit_count(o::HermitianObservable) = convert(Int, log2(size(o.matrix, 1)))
LinearAlgebra.eigvals(o::HermitianObservable) = eigvals(Hermitian(o.matrix))
unscaled(o::HermitianObservable) = o
Base.:(*)(o::HermitianObservable, n::Real) = HermitianObservable(Float64(n) .* o.matrix)
Base.show(io::IO, ho::HermitianObservable) = print(io, "HermitianObservable($(size(ho.matrix)))")

"""
    TensorProduct <: Observable
    TensorProduct(factors::Vector{<:Observable}) -> TensorProduct
    TensorProduct(factors::Vector{String}) -> TensorProduct

Struct representing a tensor product of smaller observables.

# Examples
```jldoctest
julia> BraketSimulator.Observables.TensorProduct(["x", "h"])
X @ H

julia> ho = BraketSimulator.Observables.HermitianObservable([0 1; 1 0]);

julia> BraketSimulator.Observables.TensorProduct([ho, BraketSimulator.Observables.Z()])
HermitianObservable((2, 2)) @ Z
```
"""
struct TensorProduct{O} <: Observable where {O<:Observable}
    factors::Vector{O}
    coefficient::Float64
    function TensorProduct{O}(v::Vector{O}, coefficient::Float64=1.0) where {O<:NonCompositeObservable}
        coeff = coefficient
        flattened_v = Vector{O}(undef, length(v))
        for (oi, o) in enumerate(v)
            flattened_v[oi] = unscaled(o)
            coeff *= coef(o)
        end
        return new(flattened_v, coeff)
    end
    function TensorProduct{O}(v::Vector{O}, coefficient::Float64=1.0) where {O<:Observable}
        any(v_->v_ isa Sum, v) && throw(ArgumentError("Sum observable not allowed in TensorProduct."))
        coeff = prod(coef, v, init=coefficient)
        flattened_v = Iterators.flatmap(v) do o 
            return o isa TensorProduct ? (unscaled(o) for o in o.factors) : (unscaled(o),)
        end
        flat_v = collect(flattened_v)
        return new(flat_v, coeff)
    end
end
TensorProduct(o::Vector{T}, coefficient::Float64=1.0) where {T<:Union{String, Observable}} = TensorProduct{T}(o, coefficient)
TensorProduct(o::Vector{String}, coefficient::Float64=1.0) = TensorProduct([StructTypes.constructfrom(Observable, s) for s in o], coefficient)

Base.:(*)(o::TensorProduct, n::Real) = TensorProduct(deepcopy(o.factors), Float64(n*o.coefficient))
unscaled(o::TensorProduct) = TensorProduct(o.factors, 1.0)
qubit_count(o::TensorProduct) = sum(qubit_count.(o.factors))
StructTypes.lower(x::TensorProduct{O}) where {O<:Observable} = Union{String, Vector{Vector{Vector{Float64}}}}[convert(Union{String, Vector{Vector{Vector{Float64}}}}, o) for o in mapreduce(StructTypes.lower, vcat, x.factors)]
Base.:(==)(t1::TensorProduct, t2::TensorProduct) = t1.factors == t2.factors && t1.coefficient ≈ t2.coefficient
Base.copy(t::TensorProduct) = TensorProduct(deepcopy(t.factors), t.coefficient)
_evs(os::Vector{O}, coeff::Float64=1.0) where {O<:StandardObservable} = PauliEigenvalues(Val(length(os)), coeff)
_evs(os::Vector{O}, coeff::Float64=1.0) where {O<:Observable} = mapfoldl(eigvals, kron, os, init=[coeff])::Vector{Float64}
LinearAlgebra.eigvals(o::TensorProduct{O}) where {O} = return _evs(o.factors, o.coefficient)
function Base.show(io::IO, o::TensorProduct)
    coef_str = isone(o.coefficient) ? "" : string(o.coefficient) * " * "
    print(io, coef_str)
    for f in o.factors[1:end-1]
        print(io, f)
        print(io, " @ ")
    end
    print(io, o.factors[end])
    return
end


# exclude Sum from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
"""
    Sum <: Observable
    Sum(summands::Vector{<:Observable}) -> Sum 

Struct representing the sum of observables.

# Examples
```jldoctest
julia> o1 = 2.0 * BraketSimulator.Observables.I() * BraketSimulator.Observables.Z();

julia> o2 = 3.0 * BraketSimulator.Observables.X() * BraketSimulator.Observables.X();

julia> o = o1 + o2
Sum(2.0 * I @ Z, 3.0 * X @ X)
```
"""
struct Sum <: Observable
    summands::Vector{Observable}
    coefficient::Float64
    function Sum(v)
        flattened_v = Iterators.flatmap(obs->obs isa Sum ? obs.summands : (obs,), v)
        new(collect(flattened_v), 1.0)
    end
end
Base.length(s::Sum) = length(s.summands)
Base.:(*)(s::Sum, n::Real) = Sum(Float64(n) .* deepcopy(s.summands))
function Base.:(==)(s1::Sum, s2::Sum)
    length(s1) != length(s2) && return false
    are_eq = true
    for (summand1, summand2) in zip(s1.summands, s2.summands)
        are_eq &= (summand1 == summand2)
        are_eq || return false
    end
    return true
end
function Base.:(==)(s::Sum, o::Observable)
    length(s) == 1 || return false
    return first(s.summands) == o
end
Base.:(==)(o::Observable, s::Sum) = s == o
function Base.show(io::IO, s::Sum)
    print(io, "Sum(")
    for summand in s.summands[1:end-1]
        print(io, summand)
        print(io, ", ")
    end
    print(io, s.summands[end])
    print(io, ")")
    return
end
StructTypes.lower(s::Sum) = [StructTypes.lower(summand) for summand in s.summands]

Base.:(*)(o1::O1, o2::O2) where {O1<:Observable, O2<:Observable} = TensorProduct([o1, o2], 1.0)
Base.:(*)(n::Real, o::Observable) = o*n 
Base.:(+)(o1::Observable, o2::Observable) = Sum([o1, o2])
Base.:(-)(o1::Observable, o2::Observable) = Sum([o1, -1.0 * o2])
# COV_EXCL_STOP

# nosemgrep
function StructTypes.constructfrom(::Type{Observable}, obj::String)
    (obj == "i" || obj == "I") && return I()
    (obj == "x" || obj == "X") && return X()
    (obj == "y" || obj == "Y") && return Y()
    (obj == "z" || obj == "Z") && return Z()
    (obj == "h" || obj == "H") && return H()
    throw(ArgumentError("Observable of type \"$obj\" provided, only \"i\", \"x\", \"y\", \"z\", and \"h\" are valid."))
end
StructTypes.constructfrom(::Type{Observable}, obj::Matrix{ComplexF64})              = HermitianObservable(obj)
StructTypes.constructfrom(::Type{Observable}, obj::Vector{Vector{Vector{Float64}}}) = HermitianObservable(obj)
StructTypes.constructfrom(::Type{Observable}, obj::Vector{T}) where {T} = length(obj) == 1 ? StructTypes.constructfrom(Observable, obj[1]) : TensorProduct([StructTypes.constructfrom(Observable, o) for o in obj])

end

diagonalizing_gates(g::Observables.I, targets) = Instruction[]
diagonalizing_gates(g::Observables.H, targets) =
    [Instruction(Ry(-π / 4.0), t) for t in targets]
diagonalizing_gates(g::Observables.X, targets) =
    [Instruction(H(), t) for t in targets]
diagonalizing_gates(g::Observables.Y, targets) =
    [Instruction(Unitary(1 / √2 * [1.0 -im; 1.0 im]), t) for t in targets]
diagonalizing_gates(g::Observables.Z, targets) = Instruction[]
function diagonalizing_gates(g::Observables.HermitianObservable, targets)
    size(g.matrix, 1) == 2^length(targets) &&
        return [Instruction(Unitary(eigvecs(g.matrix)), targets)]
    size(g.matrix, 1) == 2 &&
        length(targets) > 1 &&
        return [
            Instruction(Unitary(eigvecs(g.matrix)), target) for target in targets
        ]
end
diagonalizing_gates(g::Observables.TensorProduct, targets) = reduce(
    vcat,
    [diagonalizing_gates(f, t) for (f, t) in zip(g.factors, targets)],
    init = Instruction[],
)
