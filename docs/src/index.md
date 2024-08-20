```@meta
DocTestSetup = quote using BraketSimulator, BraketSimulator.Observables; using BraketSimulator: Program, Circuit, qubits, CNot, H, Rx, FreeParameter, QubitSet, AdjointGradient, BitFlip, qubit_count, Qubit, StateVector, Measure, Probability, Ry, Amplitude, Instruction, DensityMatrix, add_instruction! end
CurrentModule = BraketSimulator
```

# BraketSimulator

This package is a suite of Julia simulators of gate-based quantum circuits with (density matrix) and without (state vector) noise.
It is designed to integrate with [Amazon Braket](https://aws.amazon.com/braket/), the quantum computing service from AWS.
By default, it offers threaded CPU-based simulation of these circuits, and an optional package extension you can integrate with Python.
To use the Python integration, you will need to install [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl).

See the [Julia Pkg docs](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions)) for more information about package extensions.

If you wish to use this package **from Python**, see [`amazon-braket-simulator-v2`](), a Python package built on top of [`juliacall`](https://pypi.org/project/juliacall/)
which will automatically install Julia and all necessary Julia packages in a Python virtual environment, set appropriate environment variables, and allow you to use
these simulators from Python packages such as the [Amazon Braket SDK](https://github.com/amazon-braket/amazon-braket-sdk-python) or [PennyLane](https://pennylane.ai/).

In order to achieve the best performance for your simulations, you should set `-t auto` when you launch Julia or set the environment variable
[`JULIA_NUM_THREADS`](https://docs.julialang.org/en/v1/manual/environment-variables/#JULIA_NUM_THREADS) to `auto` (the number of CPU threads).

## Quick Start

You can install this package and its dependencies using the Julia package manager. Note that the **minimum** supported Julia version is **1.9**. If you need to install Julia itself,
follow the directions on the [JuliaLang website](https://julialang.org/downloads/).
```julia
# install the package
using Pkg
Pkg.add("BraketSimulator")
```

Then you can run a simulation of a simple [GHZ state](https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state) preparation circuit.

```jldoctest
julia> using BraketSimulator

julia> using BraketSimulator: Circuit, H, CNot, Amplitude

julia> n_qubits = 10;

julia> c = Circuit();

julia> add_instruction!(c, Instruction(H(), 0));

julia> foreach(q->add_instruction!(c, Instruction(CNot(), [0, q])), 1:n_qubits-1);

julia> push!(c.result_types, Amplitude([repeat("0", n_qubits), repeat("1", n_qubits)]));

julia> sim = StateVectorSimulator(n_qubits, 0); # use the state vector simulator (without noise)

julia> res = simulate(sim, Program(c), 0);

julia> res.resultTypes[1].value
Dict{String, ComplexF64} with 2 entries:
  "0000000000" => 0.707107+0.0im
  "1111111111" => 0.707107+0.0im
```
 
