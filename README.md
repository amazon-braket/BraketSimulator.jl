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

## Package Extensions

`BraketSimulator.jl` has optional extensions to support integration with Python (`BraketSimulatorPythonExt`) and a simulation backend based on NVIDIA's CUQUANTUM state vector simulator, `cuStateVec` (`BraketSimulatorCuStateVecExt`). The two extensions are compatible with one another.

To use `BraketSimulatorPythonExt`, you will need to install [`PythonCall.jl`](https://github.com/JuliaPy/PythonCall.jl), and then load it before `BraketSimulator.jl` like so:
```julia
using PythonCall, BraketSimulator
```
`PythonCall.jl` will try to install the necessary **Python** dependencies using the `CondaPkg.toml` present in top-level package folder. If you already have all the necessary Python dependencies installed, you can set the environment variable `JULIA_CONDAPKG_BACKEND="Null"` to have `CondaPkg.jl` use your system Python and its installed packages.

To use `BraketSimulatorCuStateVecExt`, you will need to install `cuStateVec.jl`, which is a subpackage of [`CUDA.jl`](https://github.com/JuliaGPU/CUDA.jl), and then load it before `BraketSimulator.jl` like so:
```julia
using cuStateVec, BraketSimulator
```
Currently, only NVIDIA GPUs are supported. To use the GPU-accelerated simulation backends, create a `LocalSimulator` with the backend `"braket_sv_custatevec"` (pure state simulation) or `"braket_dm_custatevec"` (noisy density matrix simulation).

## Usage Notes

If you want to use this package **from Python**, see the Python wrapper package [`amazon-braket-julia-simulator`](), which provides a Python interaction layer and will install Julia, this package and other Julia dependencies in a sandboxed virtual environment for you.

Launch Julia with the command line option `-t auto` to have Julia launch and manage threads to parallelize the simulation(s). `BraketSimulator.jl` will parallelize both across tasks within a batch and inside a single task's evolution.

Keep in mind that the first qubit has index `0`, **not** index `1`.

## Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

## License

This project is licensed under the Apache-2.0 License.
