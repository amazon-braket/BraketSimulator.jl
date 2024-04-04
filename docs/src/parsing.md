```@meta
CurrentModule = BraketSimulator
```

# OpenQASM3 Parsing

`BraketSimulators.jl` contains utilities for walking an OpenQASM3 abstract syntaxt tree (AST) based off [`OpenQASM3.jl`](). Most users will not need to interact with these, but plugin developers may wish to if they are emitting OpenQASM3 or adding support for new pragmas.

```@docs
BraketSimulator.QASMGateDefContext
BraketSimulator.lookup_def
BraketSimulator.QASMGlobalContext
BraketSimulator.QASMSubroutineContext
BraketSimulator.QASMBlockContext
```
