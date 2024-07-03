```@meta
CurrentModule = BraketSimulator
```

# Noises

`BraketSimulators.jl` provides many pre-implemented noise channels which can be applied to circuits.
Noisy circuits can be simulated used the density matrix simulator.

```@docs
BraketSimulator.Noise
BraketSimulator.Kraus
BraketSimulator.BitFlip
BraketSimulator.PhaseFlip
BraketSimulator.PauliChannel
BraketSimulator.TwoQubitPauliChannel
BraketSimulator.MultiQubitPauliChannel
BraketSimulator.Depolarizing
BraketSimulator.PhaseDamping
BraketSimulator.AmplitudeDamping
BraketSimulator.GeneralizedAmplitudeDamping
BraketSimulator.TwoQubitDepolarizing
BraketSimulator.TwoQubitDephasing
```
