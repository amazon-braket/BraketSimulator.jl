using Test, Aqua, BraketSimulator

Aqua.test_all(BraketSimulator, ambiguities=false, piracies=false, persistent_tasks = false)

@testset "BraketSimulator" begin
    for test in (
        "python_ext",
        "openqasm",
        "sv_simulator",
        "dm_simulator",
        "utils",
        "result_types",
        "braket_integration",
    )
        @testset "$test" begin
            include(test * ".jl")
        end
    end
end
