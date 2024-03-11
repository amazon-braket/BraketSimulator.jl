using BraketCircuitSimulator
using Test
using Aqua

@testset "BraketCircuitSimulator.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(BraketCircuitSimulator)
    end
    # Write your tests here.
end
