```@meta
CurrentModule = BraketSimulator
```

# Internals

These functions are for internal use for preparing IRs for simulation and validating properties to make sure all instructions and results are supported.

```@docs
BraketSimulator._combine_operations
BraketSimulator._prepare_program
BraketSimulator._get_measured_qubits
BraketSimulator._compute_results
BraketSimulator.flip_bit
BraketSimulator.flip_bits
BraketSimulator.pad_bit
BraketSimulator.pad_bits
BraketSimulator.matrix_rep
BraketSimulator.endian_qubits
BraketSimulator.get_amps_and_qubits
```
