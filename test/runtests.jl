using Test, Aqua, Documenter, BraketSimulator

Aqua.test_all(BraketSimulator, ambiguities=false)
Aqua.test_ambiguities(BraketSimulator)
dir_list = filter(x-> startswith(x, "test_") && endswith(x, ".jl"), readdir(@__DIR__))

@testset "BraketSimulator" begin
    @testset "$test" for test in dir_list
        include(test)
    end
    @testset "docs" begin
          Documenter.DocMeta.setdocmeta!(BraketSimulator, :DocTestSetup, :(using BraketSimulator, BraketSimulator.Observables; using BraketSimulator: Program, Circuit, qubits, CNot, H, Rx, FreeParameter, QubitSet, AdjointGradient, BitFlip, qubit_count, Qubit, StateVector, Measure, Probability, Ry, Amplitude, Instruction, DensityMatrix, add_instruction!); recursive=true)
          Documenter.doctest(BraketSimulator)
    end
end
