using BraketSimulator
using Documenter

DocMeta.setdocmeta!(BraketSimulator, :DocTestSetup, :(using BraketSimulator, BraketSimulator.Observables; using BraketSimulator: Program, Circuit, qubits, CNot, H, Rx, FreeParameter, QubitSet, AdjointGradient, BitFlip, qubit_count, Qubit, StateVector, Measure, Probability, Ry, Amplitude, Instruction, DensityMatrix, add_instruction!); recursive=true)

makedocs(;
    modules=[BraketSimulator],
    sitename="BraketSimulator.jl",
    repo="github.com/amazon-braket/BraketSimulator.jl",
    format=Documenter.HTML(;
        canonical="github.com/amazon-braket/BraketSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Simulators" => "sims.md",
        "Gates" => "gates.md",
        "Custom Gates" => "custom_gates.md",
        "Noises" => "noises.md",
        "Observables" => "observables.md",
        "Results" => "results.md",
    ],
)

deploydocs(;
    repo="github.com/amazon-braket/BraketSimulator.jl",
    devbranch="main",
)
