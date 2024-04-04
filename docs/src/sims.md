```@meta
CurrentModule = BraketSimulator
```

# Simulators

`BraketSimulators.jl` provides two types of simulators: `StateVectorSimulator` for pure state simulation (without noise) and `DensityMatrixSimulator` for noisy simulation.
Each type is [parameterized](https://docs.julialang.org/en/v1/manual/types/#Parametric-Types) by an **element type** (which should be a Julia `Complex` type, such as `ComplexF64`)
and an **array type** (so that we can specialize for GPU arrays, for example).

Each simulator can be initialized with a `qubit_count` and `shots` value. You may query the [`properties`](@ref Braket.properties) of a simulator to learn what gate types, result types, and other operations it supports.

```@docs
AbstractSimulator
StateVectorSimulator
DensityMatrixSimulator
Braket.properties(::StateVectorSimulator)
Braket.properties(::DensityMatrixSimulator)
BraketSimulator.evolve!
expectation
probabilities
```
