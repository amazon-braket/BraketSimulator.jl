using Test, LinearAlgebra, Braket, Braket.Observables, BraketSimulator
import Braket: Instruction, PauliEigenvalues

herm_mat = [-0.32505758-0.32505758im 0.88807383; 0.62796303+0.62796303im 0.45970084]

single_qubit_tests = [
                      (Observables.H(), (Ry(-π / 4),), PauliEigenvalues(Val(1))),
                      (Observables.X(), (H(),), PauliEigenvalues(Val(1))),
                      (Observables.Y(), (Z(), S(), H()), PauliEigenvalues(Val(1))),
                      (Observables.Z(), (), PauliEigenvalues(Val(1))),
                      (Observables.I(), (), [1, 1]),
                      (Observables.HermitianObservable([1 1-im; 1+im -1]), (Unitary(herm_mat),), [-√3, √3]),
                     ]

@testset "Observables" begin
    @testset "Single qubit $obs" for (obs, expected_gates, eigenvalues) in
                                     single_qubit_tests
        actual_gates = Braket.basis_rotation_gates(obs)
        for (actual_gate, expected_gate) in zip(actual_gates, expected_gates)
            @test BraketSimulator.matrix_rep(actual_gate) ≈ BraketSimulator.matrix_rep(expected_gate)
        end
        @test collect(eigvals(obs)) ≈ collect(eigenvalues)
    end
    @testset "Tensor product of Pauli-like (with eigenvalues ±1) observables" begin
        tensor = Observables.TensorProduct([
            Observables.H(),
            Observables.X(),
            Observables.Z(),
            Observables.Y(),
        ])
        @test eigvals(tensor) == PauliEigenvalues(Val(4))

        actual_gates = Braket.basis_rotation_gates(tensor)
        @test length(actual_gates) == 4
        @test actual_gates[1] == Braket.basis_rotation_gates(Observables.H())
        @test actual_gates[2] == Braket.basis_rotation_gates(Observables.X())
        @test actual_gates[4] == Braket.basis_rotation_gates(Observables.Y())
    end
    @testset "Tensor product of arbitrary observables" begin
        tensor = Observables.TensorProduct([
            Observables.H(),
            Observables.I(),
            Observables.X(),
            Observables.Z(),
            Observables.Y(),
        ])
        @test eigvals(tensor) == mapreduce(eigvals, kron, tensor.factors)
        actual_gates = Braket.basis_rotation_gates(tensor)
        @test length(actual_gates) == 5
        @test actual_gates[1] == Braket.basis_rotation_gates(Observables.H())
        @test actual_gates[3] == Braket.basis_rotation_gates(Observables.X())
        @test actual_gates[5] == Braket.basis_rotation_gates(Observables.Y())
    end
end
