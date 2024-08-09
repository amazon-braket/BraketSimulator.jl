```@meta
CurrentModule = BraketSimulator
```

# Gates

`BraketSimulators.jl` provides many pre-implemented gates which can be used to build up circuits. For gates with angle parameters, you can supply `Irrational`s like `π` as arguments.

```@docs
BraketSimulator.Gate
BraketSimulator.AngledGate
BraketSimulator.I
BraketSimulator.X
BraketSimulator.Y
BraketSimulator.Z
BraketSimulator.H
BraketSimulator.Rx
BraketSimulator.Ry
BraketSimulator.Rz
BraketSimulator.V
BraketSimulator.Vi
BraketSimulator.T
BraketSimulator.Ti
BraketSimulator.S
BraketSimulator.Si
BraketSimulator.U
BraketSimulator.Unitary
BraketSimulator.PhaseShift
BraketSimulator.MultiQubitPhaseShift
BraketSimulator.PRx
BraketSimulator.GPi
BraketSimulator.GPi2
BraketSimulator.XX
BraketSimulator.XY
BraketSimulator.YY
BraketSimulator.ZZ
BraketSimulator.ECR
BraketSimulator.MS
BraketSimulator.CPhaseShift
BraketSimulator.CPhaseShift00
BraketSimulator.CPhaseShift01
BraketSimulator.CPhaseShift10
BraketSimulator.CNot
BraketSimulator.CY
BraketSimulator.CZ
BraketSimulator.CV
BraketSimulator.Swap
BraketSimulator.PSwap
BraketSimulator.ISwap
BraketSimulator.CCNot
BraketSimulator.CSwap
```
