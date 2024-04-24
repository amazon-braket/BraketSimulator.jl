using Test, Aqua, BraketSimulator

Aqua.test_all(BraketSimulator, ambiguities=false, piracies=false, persistent_tasks = false)

dir_list = filter(x-> startswith(x, "test_") && endswith(x, ".jl"), readdir(@__DIR__))

@testset "BraketSimulator" begin
    for test in dir_list
        if startswith(test, "test_")
            @testset "$test" begin
                include(test)
            end
        end
    end
end
