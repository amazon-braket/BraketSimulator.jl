"""
    Result

Abstract type representing a measurement to perform on
a [`Circuit`](@ref).

See also: [`Expectation`](@ref), [`Variance`](@ref),
[`Sample`](@ref), [`Probability`](@ref),
[`DensityMatrix`](@ref), and [`Amplitude`](@ref).
"""
abstract type Result end

for (typ, ir_typ, label) in ((:Expectation, :(IR.Expectation), "expectation"), (:Variance, :(IR.Variance), "variance"), (:Sample, :(IR.Sample), "sample"))
    @eval begin
        @doc """
            $($typ) <: Result
        
        Struct which represents a $($label) measurement on a [`Circuit`](@ref). 
        """
        struct $typ <: Result
            observable::Observables.Observable
            targets::QubitSet
            @doc """
                $($typ)(o, targets) -> $($typ)
                $($typ)(o) -> $($typ)
            
            Constructs a $($typ) of an observable `o` on qubits `targets`.
            
            `o` may be one of:
              - Any [`Observable`](@ref Observables.Observable)
              - A `String` corresponding to an `Observable` (e.g. `\"x\"``)
              - A `Vector{String}` in which each element corresponds to an `Observable`

            `targets` may be one of:
              - A [`QubitSet`](@ref)
              - A `Vector` of `Int`s and/or [`Qubit`](@ref)s
              - An `Int` or `Qubit`
              - Absent, in which case the observable `o` will be applied to all qubits provided it is a single qubit observable.
            """ $typ(o::Observables.Observable, targets) = new(o, QubitSet(targets))
        end
        label(::$typ) = $label
        StructTypes.lower(x::$typ) = $ir_typ(StructTypes.lower(x.observable), (isempty(x.targets) ? nothing : Int.(x.targets)), $label)
    end
end


for (typ, ir_typ, label) in ((:Probability, :(IR.Probability), "probability"), (:DensityMatrix, :(IR.DensityMatrix), "densitymatrix"))
    @eval begin
        @doc """
            $($typ) <: Result
        
        Struct which represents a $($label) measurement on a [`Circuit`](@ref). 
        """
        struct $typ <: Result
            targets::QubitSet
            $typ(targets::QubitSet) = new(targets)
        end
        Base.:(==)(p1::$typ, p2::$typ) = (p1.targets == p2.targets)
        $typ()  = $typ(QubitSet())
        @doc """
            $($typ)(targets) -> $($typ)
            $($typ)() -> $($typ)

        Constructs a $($typ) on qubits `targets`.

        `targets` may be one of:
          - A [`QubitSet`](@ref)
          - A `Vector` of `Int`s and/or [`Qubit`](@ref)s
          - An `Int` or `Qubit`
          - Absent, in which case the measurement will be applied to all qubits.
        """ $typ(targets) = $typ(QubitSet(targets))
        $typ(targets::Vararg{IntOrQubit}) = $typ(QubitSet(targets...))
        label(::$typ) = $label
        StructTypes.lower(x::$typ) = $ir_typ(isempty(x.targets) ? nothing : Int.(x.targets), $label)
    end
end

# exclude adjoint gradient from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
"""
    AdjointGradient <: Result

Struct which represents a gradient computation using the adjoint differentiation method on a [`Circuit`](@ref). 
"""
struct AdjointGradient <: Result
    observable::Observable
    targets::Vector{QubitSet}
    parameters::Vector{String}
    function AdjointGradient(observable::Observable, targets::Vector{QubitSet}, parameters::Vector{String}=["all"])
        if observable isa Sum
            length(targets) == length(observable) || throw(DimensionMismatch("length of targets ($(length(targets))) must be the same as number of summands ($(length(observable)))."))
            all(length(term_target) == qubit_count(summand) for (term_target, summand) in zip(targets, observable.summands)) || throw(DimensionMismatch("each target must be the same size as the qubit count of its corresponding term."))
        else
            (length(targets) == 1 && length(targets[1]) == qubit_count(observable)) || throw(DimensionMismatch("targets $targets must have only one element if adjoint gradient observable is not a Sum."))
        end
        new(observable, targets, parameters)
    end
end

"""
    AdjointGradient(o::Observable, targets, parameters::Vector) -> AdjointGradient
    AdjointGradient(o::Observable, targets) -> AdjointGradient

Constructs an `AdjointGradient` with respect to the expectation value of
an observable `o` on qubits `targets`. The gradient will be calculated by
computing partial derivatives with respect to `parameters`. If `parameters`
is not present, is empty, or is `["all"]`, all parameters in the circuit
will be used. 

`targets` may be one of:
  - A [`QubitSet`](@ref)
  - A `Vector` of `Int`s and/or [`Qubit`](@ref)s
  - An `Int` or `Qubit`

`AdjointGradient` supports using [`Sum`](@ref) observables. If `o` is a `Sum`,
`targets` should be a nested vector of target qubits, such that the `n`-th term of
`targets` has the same length as the `n`-th term of `o`.

# Examples
```jldoctest
julia> α = FreeParameter(:alpha);

julia> op = 2.0 * Observables.X() * Observables.X();

julia> AdjointGradient(op, [QubitSet(0, 1)], [α]);
```

Using a `Sum`:

```jldoctest
julia> α = FreeParameter(:alpha);

julia> op1 = 2.0 * Observables.X() * Observables.X();

julia> op2 = -3.0 * Observables.Y() * Observables.Y();

julia> AdjointGradient(op1 + op2, [QubitSet(0, 1), QubitSet(0, 2)], [α]);
```
"""
function AdjointGradient(observable::Observable, targets::Vector{QubitSet}, parameters::Vector)
    isempty(parameters) && return AdjointGradient(observable, targets, ["all"])
    return AdjointGradient(observable, targets, string.(parameters))
end
AdjointGradient(observable::Observable, targets::QubitSet, parameters) = AdjointGradient(observable, [targets], parameters)
AdjointGradient(observable::Observable, targets::Vector{Vector{T}}, args...) where {T} = AdjointGradient(observable, [QubitSet(t) for t in targets], args...)
AdjointGradient(observable::Observable, targets::Vector{<:IntOrQubit}, args...) = AdjointGradient(observable, [QubitSet(targets)], args...)
AdjointGradient(observable::Observable, targets::IntOrQubit, args...) = AdjointGradient(observable, [QubitSet(targets)], args...)

StructTypes.StructType(::Type{AdjointGradient}) = StructTypes.CustomStruct()
function StructTypes.lower(x::AdjointGradient)
    lowered_obs     = StructTypes.lower(x.observable)
    lowered_targets = (isempty(x.targets) ? nothing : convert(Vector{Vector{Int}}, x.targets))
    IR.AdjointGradient(x.parameters, lowered_obs, lowered_targets, "adjoint_gradient")
end
# COV_EXCL_STOP

"""
    Amplitude <: Result

Struct which represents an amplitude measurement on a [`Circuit`](@ref). 
"""
struct Amplitude <: Result
    states::Vector{String}
end
"""
    Amplitude(states) -> Amplitude

Constructs an Amplitude measurement of `states`.

`states` may be one of:
  - A `Vector{String}`
  - A `String`
All elements of `states` must be `'0'` or `'1'`.
"""
Amplitude(s::String) = Amplitude([s])
Base.:(==)(a1::Amplitude, a2::Amplitude) = (a1.states == a2.states)
label(::Amplitude) = "amplitude"

"""
    StateVector <: Result

Struct which represents a state vector measurement on a [`Circuit`](@ref). 
"""
struct StateVector <: Result end
Base.:(==)(sv1::StateVector, sv2::StateVector) = true
label(::StateVector) = "statevector"

const ObservableResult = Union{Expectation, Variance, Sample}
const ObservableParameterResult = Union{AdjointGradient,}
StructTypes.lower(x::Amplitude) = IR.Amplitude(x.states, "amplitude")
StructTypes.lower(x::StateVector) = IR.StateVector("statevector")
