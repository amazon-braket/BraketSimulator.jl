```@meta
CurrentModule = BraketSimulator
```

# BraketSimulator

Documentation for [BraketSimulator](https://github.com/amazon-braket/BraketSimulator.jl).

```@index
StateVectorSimulator{T, S<:AbstractVector{T}}
StateVectorSimulator([::T], qubit_count::Int, shots::Int)
properties(svs::StateVectorSimulator) -> GateModelSimulatorDeviceCapabilities
evolve!(svs::StateVectorSimulator{T, S<:AbstractVector{T}}, operations::Vector{Instruction})
expectation
probabilities
```

```@autodocs
Modules = [BraketSimulator]
```
