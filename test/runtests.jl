using Test, Aqua, BraketSimulator

Aqua.test_all(BraketSimulator, ambiguities=false)
Aqua.test_ambiguities(BraketSimulator)
dir_list = filter(x-> startswith(x, "test_") && endswith(x, ".jl"), readdir(@__DIR__))

@testset "BraketSimulator" begin
    @testset "$test" for test in dir_list
        include(test)
    end
end
