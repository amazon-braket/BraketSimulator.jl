using BraketSimulator
using Documenter

DocMeta.setdocmeta!(BraketSimulator, :DocTestSetup, :(using BraketSimulator); recursive=true)

makedocs(;
    modules=[BraketSimulator],
    authors="Katharine Hyatt <hyatkath@amazon.com> and contributors",
    sitename="BraketSimulator.jl",
    format=Documenter.HTML(;
        canonical="https://amazon-braket.github.io/BraketSimulator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/amazon-braket/BraketSimulator.jl",
    devbranch="main",
)
