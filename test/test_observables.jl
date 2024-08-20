using Test, LinearAlgebra, BraketSimulator.Observables, BraketSimulator

herm_mat = [-0.32505758-0.32505758im 0.88807383; 0.62796303+0.62796303im 0.45970084]

single_qubit_tests = [
                      (Observables.H(), (BraketSimulator.Ry(-π / 4),), BraketSimulator.PauliEigenvalues(Val(1))),
                      (Observables.X(), (BraketSimulator.H(),), BraketSimulator.PauliEigenvalues(Val(1))),
                      (Observables.Y(), (BraketSimulator.Z(), BraketSimulator.S(), BraketSimulator.H()), BraketSimulator.PauliEigenvalues(Val(1))),
                      (Observables.Z(), (), BraketSimulator.PauliEigenvalues(Val(1))),
                      (Observables.I(), (), [1, 1]),
                      (Observables.HermitianObservable([1.0 1.0-im; 1.0+im -1.0]), (BraketSimulator.Unitary(herm_mat),), [-√3, √3]),
                     ]

@testset "Observables" begin
    @testset "Single qubit $obs" for (obs, expected_gates, eigenvalues) in
                                     single_qubit_tests
        actual_gates = BraketSimulator.basis_rotation_gates(obs)
        for (actual_gate, expected_gate) in zip(actual_gates, expected_gates)
            @test BraketSimulator.matrix_rep(actual_gate) ≈ BraketSimulator.matrix_rep(expected_gate)
        end
        @test collect(eigvals(obs)) ≈ collect(eigenvalues)
        @test BraketSimulator.qubit_count(obs) == 1
        @test ishermitian(obs)
        @test copy(obs) == obs
        if !(obs isa Observables.HermitianObservable)
            @test BraketSimulator.Observables.unscaled(2.0 * obs) == obs
        else
            scaled = 2.0 * obs
            @test scaled.matrix == obs.matrix .* 2.0
        end
    end
    @test_throws ArgumentError("Observable of type \"c\" provided, only \"i\", \"x\", \"y\", \"z\", and \"h\" are valid.") BraketSimulator.StructTypes.constructfrom(Observables.Observable, "c")
    @testset "Tensor product of BraketSimulator.Pauli-like (with eigenvalues ±1) observables" begin
        tensor = Observables.TensorProduct(["h", "x", "z", "y"])
        @test eigvals(tensor) == BraketSimulator.PauliEigenvalues(Val(4))
        @test BraketSimulator.StructTypes.constructfrom(Observables.Observable, ["h", "x", "z", "y"]) == tensor
        @test eigvals(tensor)[[1, 2]] == [1, -1]

        actual_gates = BraketSimulator.basis_rotation_gates(tensor)
        @test length(actual_gates) == 4
        @test actual_gates[1] == BraketSimulator.basis_rotation_gates(Observables.H())
        @test actual_gates[2] == BraketSimulator.basis_rotation_gates(Observables.X())
        @test actual_gates[4] == BraketSimulator.basis_rotation_gates(Observables.Y())
        @test BraketSimulator.qubit_count(tensor) == 4
        @test ishermitian(tensor)
        @test copy(tensor) == tensor
        @test BraketSimulator.Observables.unscaled(2.0 * tensor) == tensor
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
        actual_gates = BraketSimulator.basis_rotation_gates(tensor)
        @test length(actual_gates) == 5
        @test actual_gates[1] == BraketSimulator.basis_rotation_gates(Observables.H())
        @test actual_gates[3] == BraketSimulator.basis_rotation_gates(Observables.X())
        @test actual_gates[5] == BraketSimulator.basis_rotation_gates(Observables.Y())
        @test ishermitian(tensor)
        @test copy(tensor) == tensor
        @test BraketSimulator.Observables.unscaled(2.0 * tensor) == tensor
    end
end
