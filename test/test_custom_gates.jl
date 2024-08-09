using Test, Logging, BraketSimulator, DataStructures

@testset "Custom gates" begin
    @testset "Double excitation" begin
        ϕ  = 3.56
        nq = 4
        @test BraketSimulator.qubit_count(BraketSimulator.DoubleExcitation(0.1)) == 4
        instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]), BraketSimulator.Instruction(BraketSimulator.DoubleExcitation(ϕ), [0, 1, 2, 3])]
        # instructions for the gate decomposition (from PennyLane)
        de_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 2]), BraketSimulator.Instruction(BraketSimulator.H(), [3]), BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]), BraketSimulator.Instruction(BraketSimulator.Ry(ϕ/8), [1]), BraketSimulator.Instruction(BraketSimulator.Ry(-ϕ/8), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 3]), BraketSimulator.Instruction(BraketSimulator.H(), [3]), BraketSimulator.Instruction(BraketSimulator.CNot(), [3, 1]), BraketSimulator.Instruction(BraketSimulator.Ry(ϕ/8), [1]), BraketSimulator.Instruction(BraketSimulator.Ry(-ϕ/8), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 1]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 0]), BraketSimulator.Instruction(BraketSimulator.Ry(-ϕ/8), [1]), BraketSimulator.Instruction(BraketSimulator.Ry(ϕ/8), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [3, 1]), BraketSimulator.Instruction(BraketSimulator.H(), [3]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 3]), BraketSimulator.Instruction(BraketSimulator.Ry(-ϕ/8), [1]), BraketSimulator.Instruction(BraketSimulator.Ry(ϕ/8), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 0]), BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [3]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 2]), BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3])]
        # instructions for the matrix representation
        u_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]), BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(BraketSimulator.DoubleExcitation(ϕ)))), [0, 1, 2, 3])]
        # initial state
        # [ 1/√2, 0, 0, 1/√2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # final state
        # [ 1/√2, 0, 0, -sin(ϕ/2)/√2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)/√2, 0, 0, 0]
        state_vector = [1/√2, 0, 0, -sin(ϕ/2)/√2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)/√2, 0, 0, 0]
        probability_amplitudes = 0.5*[1, 0, 0, sin(ϕ/2)^2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)^2, 0, 0, 0]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = BraketSimulator.evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test BraketSimulator.qubit_count(BraketSimulator.DoubleExcitation(ϕ)) == 4
        @test BraketSimulator.qubit_count(BraketSimulator.DoubleExcitation) == 4
        @test inv(BraketSimulator.DoubleExcitation(ϕ)) == BraketSimulator.DoubleExcitation(ϕ, -1.0)
    end
    @testset "Single excitation" begin
        ϕ  = 3.56
        nq = 2
        @test BraketSimulator.qubit_count(BraketSimulator.SingleExcitation(0.1)) == 2
        instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.SingleExcitation(ϕ), [0, 1])]
        de_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.Ti(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.S(), [0]), BraketSimulator.Instruction(BraketSimulator.Ti(), [1]), BraketSimulator.Instruction(BraketSimulator.Si(), [1]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.CNot(), [1, 0]), BraketSimulator.Instruction(BraketSimulator.Rz(-ϕ/2), [0]), BraketSimulator.Instruction(BraketSimulator.Ry(ϕ/2), [1]), BraketSimulator.Instruction(BraketSimulator.CNot(), [1, 0]), BraketSimulator.Instruction(BraketSimulator.Si(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.T(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.S(), [1]), BraketSimulator.Instruction(BraketSimulator.T(), [1])]
        # instructions for the matrix representation
        u_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(BraketSimulator.SingleExcitation(ϕ)))), [0, 1])]
        # initial state
        # [ 1/2, 1/2, 1/2, 1/2]
        # final state
        # [ 1/2, cos(ϕ/2) -/2, 1/2, 1/2 ]
        state_vector = 0.5 * [1, cos(ϕ/2) - sin(ϕ/2), cos(ϕ/2) + sin(ϕ/2), 1]
        probability_amplitudes = 0.25*[1, (cos(ϕ/2) - sin(ϕ/2))^2, (cos(ϕ/2) + sin(ϕ/2))^2, 1]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = BraketSimulator.evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test BraketSimulator.qubit_count(BraketSimulator.SingleExcitation(ϕ)) == 2
        @test BraketSimulator.qubit_count(BraketSimulator.SingleExcitation) == 2
        @test inv(BraketSimulator.SingleExcitation(ϕ)) == BraketSimulator.SingleExcitation(ϕ, -1.0)
    end
    @testset "3-angle U" begin
        θ = 1.34
        ϕ = 2.12
        λ = 0.43
        nq = 2
        @test BraketSimulator.qubit_count(BraketSimulator.U(θ, ϕ, λ)) == 1
        instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.U(θ, ϕ, λ), [0])]
        u_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(BraketSimulator.U(θ, ϕ, λ)))), [0])]
        # initial state
        # [ 1/2, 1/2, 1/2, 1/2]
        # final state
        # 1/2 * [ cos(θ/2) - exp(im*λ)*sin(θ/2), exp(im*ϕ)*sin(θ/2) + exp(im*(ϕ+λ))*cos(ϕ/2),  cos(θ/2) - exp(im*λ)*sin(θ/2), exp(im*ϕ)*sin(θ/2) + exp(im*(ϕ+λ))*cos(ϕ/2)]
        state_vector = ComplexF64[0.10968333631828209 - 0.12943546335661305im, 0.10968333631828209 - 0.12943546335661305im, -0.4873868533239797 + 0.483394333610458im, -0.4873868533239797 + 0.483394333610458im] 
        probability_amplitudes = [0.0287839734402505, 0.0287839734402505, 0.47121602655974926, 0.47121602655974926] 
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = BraketSimulator.evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test BraketSimulator.qubit_count(BraketSimulator.U(θ, ϕ, λ)) == 1
        @test BraketSimulator.qubit_count(BraketSimulator.U) == 1
        @test inv(BraketSimulator.U(θ, ϕ, λ)) == BraketSimulator.U(θ, ϕ, λ, -1)
    end
    @testset "MultiQubitPhaseShift" begin
        ϕ = 2.12
        @testset "Simulator $sim, n_qubits $nq" for sim in (StateVectorSimulator, DensityMatrixSimulator), nq in 1:4
            instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ), collect(0:nq-1))])
            u_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ)))), collect(0:nq-1))])
            state_vector = exp(im*ϕ)/√(2^nq) * ones(2^nq)
            probability_amplitudes = 1/2^nq * ones(2^nq)
            @test BraketSimulator.qubit_count(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ)) == nq
            @testset "BraketSimulator.Instruction set $ix_label" for (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
                simulation = sim(nq, 0)
                simulation = BraketSimulator.evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
        @test BraketSimulator.qubit_count(BraketSimulator.MultiQubitPhaseShift{2}(ϕ)) == 2
        @test BraketSimulator.qubit_count(BraketSimulator.MultiQubitPhaseShift{3}(ϕ)) == 3
        @test inv(BraketSimulator.MultiQubitPhaseShift{2}(ϕ)) == BraketSimulator.MultiQubitPhaseShift{2}(ϕ, -1.0)
    end
    @testset "MultiRz" begin
        ϕ = 2.12
        @testset "Simulator $sim, n_qubits $nq" for sim in (StateVectorSimulator, DensityMatrixSimulator), nq in 1:6
            instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.MultiRZ(ϕ), collect(0:nq-1))])
            simulation = sim(nq, 0)
            simulation = BraketSimulator.evolve!(simulation, instructions)
            state_vector = √(1/2^nq) * exp.(-im*ϕ/2 .* BraketSimulator.PauliEigenvalues(Val(nq)))
            probability_amplitudes = 1/2^nq * ones(2^nq)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test BraketSimulator.qubit_count(BraketSimulator.MultiRZ(ϕ)) == 1
        @test inv(BraketSimulator.MultiRZ(ϕ)) == BraketSimulator.MultiRZ(ϕ, -1.0)
    end
    @testset "Control" begin
        x = BraketSimulator.Control(BraketSimulator.X(), ())
        @test BraketSimulator.matrix_rep(x) == BraketSimulator.matrix_rep(BraketSimulator.X())
        @test BraketSimulator.qubit_count(x) == 1
        cx = BraketSimulator.Control(BraketSimulator.X(), (1,))
        @test Matrix(BraketSimulator.matrix_rep(cx)) == [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0] 
        @test BraketSimulator.qubit_count(cx) == 2
        cx = BraketSimulator.Control(BraketSimulator.Control(BraketSimulator.X(), (1,)), ())
        @test Matrix(BraketSimulator.matrix_rep(cx)) == [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0] 
        @test BraketSimulator.qubit_count(cx) == 2
        ccx = BraketSimulator.Control(cx, (1,)) 
        @test BraketSimulator.qubit_count(ccx) == 3
        @testset for (nq, g, cg) in ((2, cx, BraketSimulator.CNot()), (3, ccx, BraketSimulator.CCNot()))
            instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [0])], [BraketSimulator.Instruction(BraketSimulator.CNot(), [q, q+1]) for q in 0:nq-3], [BraketSimulator.Instruction(g, collect(0:nq-1))])
            canonical_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [0])], [BraketSimulator.Instruction(BraketSimulator.CNot(), [q, q+1]) for q in 0:nq-3], [BraketSimulator.Instruction(cg, collect(0:nq-1))])
            vec = zeros(2^nq)
            vec[1] = 1
            vec[end] = 1
            state_vector = 1/√2 * vec 
            probability_amplitudes = 1/2 * vec 
            @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
                (ix_label, ixs) in (("control", instructions), ("builtin", canonical_instructions))
                simulation = sim(nq, 0)
                simulation = BraketSimulator.evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
        @testset for nq in 1:2
            ϕ = 2.12
            @test BraketSimulator.qubit_count(BraketSimulator.Control(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq))) == nq
            instructions     = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq)), collect(0:nq-1))])
            u                = Matrix(BraketSimulator.matrix_rep(BraketSimulator.Control(BraketSimulator.MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq))))
            u_instructions   = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Unitary(u), collect(0:nq-1))])
            state_vector     = 1/√(2^nq) * ones(ComplexF64, 2^nq)
            state_vector[1] *= exp(im*ϕ)
            probability_amplitudes = 1/2^nq * ones(2^nq)
            @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator), (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
                simulation = sim(nq, 0)
                simulation = BraketSimulator.evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
    end
end
