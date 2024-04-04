using Braket, BraketSimulator
using Documenter

DocMeta.setdocmeta!(BraketSimulator, :DocTestSetup, :(using Braket, BraketSimulator); recursive=true)

makedocs(;
    modules=[BraketSimulator],
    authors="Katharine Hyatt <hyatkath@amazon.com> and contributors",
    sitename="BraketSimulator.jl",
    repo="git@ssh.gitlab.aws.dev:braket-science/braketsimulator.git",
    format=Documenter.HTML(;
        #canonical="https://amazon-braket.github.io/BraketSimulator.jl",
        canonical="git@ssh.gitlab.aws.dev:braket-science/braketsimulator.git",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Simulators" => "sims.md",
        "OpenQASM3 Parsing" => "parsing.md",
    ],
)

# TODO UNCOMMENT ME WHEN PUBLIC
#=deploydocs(;
    repo="github.com/amazon-braket/BraketSimulator.jl",
    devbranch="main",
)=#
