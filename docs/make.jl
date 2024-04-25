using Braket, BraketSimulator
using Documenter

DocMeta.setdocmeta!(BraketSimulator, :DocTestSetup, :(using Braket, BraketSimulator); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/amazon-braket/BraketSimulator.jl",
    devbranch="main",
)
