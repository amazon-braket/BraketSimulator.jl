# BraketSimulator

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://amazon-braket.github.io/BraketSimulator.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://amazon-braket.github.io/BraketSimulator.jl/dev/)
[![Build Status](https://github.com/amazon-braket/BraketSimulator.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/amazon-braket/BraketSimulator.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/amazon-braket/BraketSimulator.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/amazon-braket/BraketSimulator.jl)
[![Coverage](https://coveralls.io/repos/github/amazon-braket/BraketSimulator.jl/badge.svg?branch=main)](https://coveralls.io/github/amazon-braket/BraketSimulator.jl?branch=main)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package is a quantum circuit simulator written in the [Julia programming language](https://julialang.org/), meant to be compatible with the [Amazon Braket SDK](https://github.com/aws/amazon-braket-sdk-python). It can simulate gate-based quantum circuits using both statevectors and density matrices (when noise is present).

## Installation & Prerequisites

You do *not* need a Python installation or the Python Amazon Braket SDK installed to use this package.

**Supported operating systems and architectures:**

| CPU arch \ OS | Linux | macOS | Windows |
| ------------- | ----- | ----- | ------- |
| x86\_64       | :white_check_mark: | :white_check_mark: | :white_check_mark: |
| aarch64 (ARM) | :white_check_mark: | :white_check_mark: | :interrobang: |


All necessary Julia packages will be installed for you when you run `Pkg.add("BraketSimulator")` or `] instantiate` (if you're doing a `dev` install). The minimum supported Julia version is 1.9.

## Generating precompilation instructions

In order to generate a new set of precompilation instructions, you'll need to comment out the `include("precompile.jl")` lines in `src/BraketSimulator.jl` and `ext/BraketSimulatorPythonExt/BraketSimulatorPythonExt.jl`. Then, `cd` to `precompile` and run:

```bash
julia --project=../test snoop_compilation.jl
```

This will run the test suite using [`SnoopCompile.jl`](https://timholy.github.io/SnoopCompile.jl/dev/snoopi_deep_parcel/) to generate precompilation instructions. Once this is completed, you can *un*comment `include("precompile.jl")` lines in `src/BraketSimulator.jl` and `ext/BraketSimulatorPythonExt/BraketSimulatorPythonExt.jl`.

## Package Extensions

`BraketSimulator.jl` has an optional extension to support integration with Python (`BraketSimulatorPythonExt`).

To use `BraketSimulatorPythonExt`, you will need to install [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl), and then load it before `BraketSimulator.jl` like so:
```julia
using PythonCall, BraketSimulator
```
`PythonCall.jl` will try to install the necessary **Python** dependencies using the `CondaPkg.toml` present in top-level package folder. If you already have all the necessary Python dependencies installed, you can set the environment variable `JULIA_CONDAPKG_BACKEND="Null"` to have `CondaPkg.jl` use your system Python and its installed packages.

## Usage Notes

If you want to use this package **from Python**, see the Python wrapper package [`amazon-braket-julia-simulator`](), which provides a Python interaction layer and will install Julia, this package and other Julia dependencies in a sandboxed virtual environment for you.

Launch Julia with the command line option `-t auto` to have Julia launch and manage threads to parallelize the simulation(s). `BraketSimulator.jl` will parallelize both across tasks within a batch and inside a single task's evolution.

Keep in mind that the first qubit has index `0`, **not** index `1`.

## Example

In this example, we create a [Greenberger-Horne-Zeilinger (GHZ) state](https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state) and simulate it:

```julia
using Braket, BraketSimulator

function ghz_circuit(qubit_count::Int)
    ghz_circ = Circuit()
    H(ghz_circ, 0)
    for target_qubit = 1:qubit_count-1
        CNot(ghz_circ, 0, target_qubit)
    end
    return ghz_circ 
end

device = LocalSimulator("braket_sv_v2")
n_qubits = 10
task_result = simulate(device, ghz_circuit(n_qubits), shots=10)
```

## Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

## License

This project is licensed under the Apache-2.0 License.
