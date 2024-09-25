```@meta
CurrentModule = BraketSimulator
```

# Circuits

Circuits are made up of *instructions* (operations to apply to the qubits -- [gates](gates.md) and [noises](noises.md)) and *result types* ([results](results.md)).
OpenQASM3 programs are parsed to circuits which are then run on the simulator.

```@docs
BraketSimulator.Circuit
BraketSimulator.Operator
BraketSimulator.QuantumOperator
BraketSimulator.FreeParameter
BraketSimulator.Reset
BraketSimulator.Barrier
BraketSimulator.Delay
BraketSimulator.Measure
BraketSimulator.Instruction
BraketSimulator.QubitSet
BraketSimulator.qubit_count
BraketSimulator.qubits
BraketSimulator.basis_rotation_instructions!
```
