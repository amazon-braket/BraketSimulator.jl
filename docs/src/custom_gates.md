```@meta
CurrentModule = BraketSimulator
```

# Custom gates

`BraketSimulator.jl` defines some custom gates to extend what `Braket.jl` provides. These include the OpenQASM3 built-in gates `gphase` (as `MultiQubitPhaseShift`) and `U` (the single qubit three angle unitary).

```@docs
U
MultiQubitPhaseShift
MultiRZ
```
