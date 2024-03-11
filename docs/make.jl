using BraketCircuitSimulator
using Documenter

DocMeta.setdocmeta!(BraketCircuitSimulator, :DocTestSetup, :(using BraketCircuitSimulator); recursive=true)

makedocs(;
    modules=[BraketCircuitSimulator],
    authors="Katharine Hyatt <hyatkath@amazon.com> and contributors",
    sitename="BraketCircuitSimulator.jl",
    format=Documenter.HTML(;
        canonical="https://kshyatt-aws.github.io/BraketCircuitSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kshyatt-aws/BraketCircuitSimulator.jl",
    devbranch="main",
)
