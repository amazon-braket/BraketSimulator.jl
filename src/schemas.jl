"""
    Instruction
    Instruction(o::Operator, target)

Represents a single operation applied to a [`Circuit`](@ref).
Contains an `operator`, which may be any subtype of [`Operator`](@ref),
and a `target` set of qubits to which the `operator` is applied.

# Examples
```jldoctest
julia> Instruction(H(), 1)
Instruction{H}(H(1.0), QubitSet(1))

julia> Instruction(CNot(), [1, Qubit(4)])
Instruction{CNot}(CNot(1.0), QubitSet(1, 4))
```
"""
struct Instruction{O<:Operator}
    operator::O
    target::QubitSet
end
Instruction(o::O, target) where {O<:Operator} = Instruction{O}(o, QubitSet(target...))
Base.:(==)(ix1::Instruction{O}, ix2::Instruction{O}) where {O<:Operator} = (ix1.operator == ix2.operator && ix1.target == ix2.target)
bind_value!(ix::Instruction{O}, param_values::Dict{Symbol, <:Real}) where {O<:Operator} = Instruction{O}(bind_value!(ix.operator, param_values), ix.target)
remap(@nospecialize(ix::Instruction{O}), mapping::Dict{<:Integer, <:Integer}) where {O} = Instruction{O}(copy(ix.operator), [mapping[q] for q in ix.target])
