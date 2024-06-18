using Test, Logging, Braket, BraketSimulator, DataStructures

using Braket: Instruction
using BraketSimulator: DoubleExcitation, DoubleExcitationPlus, DoubleExcitationMinus, SingleExcitation, SingleExcitationPlus, SingleExcitationMinus,  matrix_rep, Control, FermionicSWAP

@testset "Custom gates" begin
    @testset "Double excitation" begin
        ϕ  = 3.56
        nq = 4
        instructions = [Instruction(H(), [0]), Instruction(CNot(), [0, 1]), Instruction(DoubleExcitation(ϕ), [0, 1, 2, 3])]
        # instructions for the gate decomposition (from PennyLane)
        de_instructions = [Instruction(H(), [0]), Instruction(CNot(), [0, 1]), Instruction(CNot(), [2, 3]), Instruction(CNot(), [0, 2]), Instruction(H(), [3]), Instruction(H(), [0]), Instruction(CNot(), [2, 3]), Instruction(CNot(), [0, 1]), Instruction(Ry(ϕ/8), [1]), Instruction(Ry(-ϕ/8), [0]), Instruction(CNot(), [0, 3]), Instruction(H(), [3]), Instruction(CNot(), [3, 1]), Instruction(Ry(ϕ/8), [1]), Instruction(Ry(-ϕ/8), [0]), Instruction(CNot(), [2, 1]), Instruction(CNot(), [2, 0]), Instruction(Ry(-ϕ/8), [1]), Instruction(Ry(ϕ/8), [0]), Instruction(CNot(), [3, 1]), Instruction(H(), [3]), Instruction(CNot(), [0, 3]), Instruction(Ry(-ϕ/8), [1]), Instruction(Ry(ϕ/8), [0]), Instruction(CNot(), [0, 1]), Instruction(CNot(), [2, 0]), Instruction(H(), [0]), Instruction(H(), [3]), Instruction(CNot(), [0, 2]), Instruction(CNot(), [2, 3])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(CNot(), [0, 1]), Instruction(Unitary(Matrix(matrix_rep(DoubleExcitation(ϕ)))), [0, 1, 2, 3])]
        # initial state
        # [ 1/√2, 0, 0, 1/√2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # final state
        # [ 1/√2, 0, 0, -sin(ϕ/2)/√2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)/√2, 0, 0, 0]
        state_vector = [1/√2, 0, 0, -sin(ϕ/2)/√2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)/√2, 0, 0, 0]
        probability_amplitudes = 0.5*[1, 0, 0, sin(ϕ/2)^2, 0, 0, 0, 0, 0, 0, 0, 0, cos(ϕ/2)^2, 0, 0, 0]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(DoubleExcitation(ϕ)) == 4
        @test inv(DoubleExcitation(ϕ)) == DoubleExcitation(-ϕ)
        @test DoubleExcitation(ϕ) ^ 0 == DoubleExcitation(0.0)
        @test DoubleExcitation(ϕ) ^ 2 == DoubleExcitation(2*ϕ)
        @test DoubleExcitation(ϕ) ^ -1 == DoubleExcitation(-ϕ)
        @test DoubleExcitation(ϕ) ^ - 3 == DoubleExcitation(-3*ϕ)
        
    end
    @testset "Single excitation" begin
        ϕ  = 3.56
        nq = 2
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(SingleExcitation(ϕ), [0, 1])]
        de_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Ti(), [0]), Instruction(H(), [0]), Instruction(S(), [0]), Instruction(Ti(), [1]), Instruction(Si(), [1]), Instruction(H(), [1]), Instruction(CNot(), [1, 0]), Instruction(Rz(-ϕ/2), [0]), Instruction(Ry(ϕ/2), [1]), Instruction(CNot(), [1, 0]), Instruction(Si(), [0]), Instruction(H(), [0]), Instruction(T(), [0]), Instruction(H(), [1]), Instruction(S(), [1]), Instruction(T(), [1])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(SingleExcitation(ϕ)))), [1, 0])]
        # initial state
        # [ 1/2, 1/2, 1/2, 1/2]
        # final state
        # [ 1/2, cos(ϕ/2) -/2, 1/2, 1/2 ]
        state_vector = 0.5 * [1, cos(ϕ/2) - sin(ϕ/2), cos(ϕ/2) + sin(ϕ/2), 1]
        probability_amplitudes = 0.25*[1, (cos(ϕ/2) - sin(ϕ/2))^2, (cos(ϕ/2) + sin(ϕ/2))^2, 1]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(SingleExcitation(ϕ)) == 2
        @test inv(SingleExcitation(ϕ)) == SingleExcitation(-ϕ)
        @test SingleExcitation(ϕ) ^ 0 == SingleExcitation(0.0)
        @test SingleExcitation(ϕ) ^ 2 == SingleExcitation(2*ϕ)
        @test SingleExcitation(ϕ) ^ -1 == SingleExcitation(-ϕ)
        @test SingleExcitation(ϕ) ^ - 3 == SingleExcitation(-3*ϕ)
    end
    @testset "3-angle U" begin
        θ = 1.34
        ϕ = 2.12
        λ = 0.43
        nq = 2
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(U(θ, ϕ, λ), [0])]
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(U(θ, ϕ, λ)))), [0])]
        # initial state
        # [ 1/2, 1/2, 1/2, 1/2]
        # final state
        # 1/2 * [ cos(θ/2) - exp(im*λ)*sin(θ/2), exp(im*ϕ)*sin(θ/2) + exp(im*(ϕ+λ))*cos(ϕ/2),  cos(θ/2) - exp(im*λ)*sin(θ/2), exp(im*ϕ)*sin(θ/2) + exp(im*(ϕ+λ))*cos(ϕ/2)]
        state_vector = ComplexF64[0.10968333631828209 - 0.12943546335661305im, 0.10968333631828209 - 0.12943546335661305im, -0.4873868533239797 + 0.483394333610458im, -0.4873868533239797 + 0.483394333610458im] 
        probability_amplitudes = [0.0287839734402505, 0.0287839734402505, 0.47121602655974926, 0.47121602655974926] 
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(U(θ, ϕ, λ)) == 1
        @test inv(U(θ, ϕ, λ)) == U(-θ, -λ, -ϕ)
    end
    @testset "MultiQubitPhaseShift" begin
        ϕ = 2.12
        @testset "Simulator $sim, n_qubits $nq" for sim in (StateVectorSimulator, DensityMatrixSimulator), nq in 1:4
            instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(MultiQubitPhaseShift{nq}(ϕ), collect(0:nq-1))])
            u_instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(Unitary(Matrix(matrix_rep(MultiQubitPhaseShift{nq}(ϕ)))), collect(0:nq-1))])
            state_vector = exp(im*ϕ)/√(2^nq) * ones(2^nq)
            probability_amplitudes = 1/2^nq * ones(2^nq)
            @testset "Instruction set $ix_label" for (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
                simulation = sim(nq, 0)
                simulation = evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
        @test qubit_count(MultiQubitPhaseShift{2}(ϕ)) == 2
        @test qubit_count(MultiQubitPhaseShift{3}(ϕ)) == 3
        @test inv(MultiQubitPhaseShift{2}(ϕ)) == MultiQubitPhaseShift{2}(-ϕ)
    end
    @testset "MultiRz" begin
        ϕ = 2.12
        @testset "Simulator $sim, n_qubits $nq" for sim in (StateVectorSimulator, DensityMatrixSimulator), nq in 1:4
            instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(MultiRZ(ϕ), collect(0:nq-1))])
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, instructions)
            state_vector = √(1/2^nq) * exp.(-im*ϕ/2 .* Braket.PauliEigenvalues(Val(nq)))
            probability_amplitudes = 1/2^nq * ones(2^nq) 
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(MultiRZ(ϕ)) == 1
        @test inv(MultiRZ(ϕ)) == MultiRZ(-ϕ)
    end
    @testset "Control" begin
        x = Control(X(), ())
        @test matrix_rep(x) == matrix_rep(X())
        @test qubit_count(x) == 1
        cx = Control(X(), (1,))
        @test Matrix(matrix_rep(cx)) == [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0] 
        @test qubit_count(cx) == 2
        ccx = Control(cx, (1,)) 
        @test qubit_count(ccx) == 3
        @testset for (nq, g, cg) in ((2, cx, CNot()), (3, ccx, CCNot()))
            instructions = vcat([Instruction(H(), [0])], [Instruction(CNot(), [q, q+1]) for q in 0:nq-3], [Instruction(g, collect(0:nq-1))])
            canonical_instructions = vcat([Instruction(H(), [0])], [Instruction(CNot(), [q, q+1]) for q in 0:nq-3], [Instruction(cg, collect(0:nq-1))])
            vec = zeros(2^nq)
            vec[1] = 1
            vec[end] = 1
            state_vector = 1/√2 * vec 
            probability_amplitudes = 1/2 * vec 
            @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
                (ix_label, ixs) in (("control", instructions), ("builtin", canonical_instructions))
                simulation = sim(nq, 0)
                simulation = evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
        @testset for nq in 1:2
            ϕ = 2.12
            @test qubit_count(Control(MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq))) == nq
            instructions     = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(Control(MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq)), collect(0:nq-1))])
            u                = Matrix(matrix_rep(Control(MultiQubitPhaseShift{nq}(ϕ), ntuple(i->0, nq))))
            u_instructions   = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(Unitary(u), collect(0:nq-1))])
            state_vector     = 1/√(2^nq) * ones(ComplexF64, 2^nq)
            state_vector[1] *= exp(im*ϕ)
            probability_amplitudes = 1/2^nq * ones(2^nq)
            @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator), (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
                simulation = sim(nq, 0)
                simulation = evolve!(simulation, ixs)
                if sim == StateVectorSimulator
                    @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
                end
                @test probability_amplitudes ≈
                      collect(BraketSimulator.probabilities(simulation))
            end
        end
    end

    @testset "Single excitation plus" begin
        ϕ  = 3.56
        nq = 2
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(SingleExcitationPlus(ϕ), [0, 1])]
        # instructions for the gate decomposition (from PennyLane)
        de_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(X(), [0]), Instruction(X(), [1]), Instruction(CPhaseShift(ϕ/2), [1, 0]), Instruction(X(), [0]), Instruction(X(), [1]), Instruction(CPhaseShift(ϕ/2), [0, 1]), Instruction(CNot(), [0, 1]), Instruction(Ry(ϕ/2), [0]), Instruction(CNot(), [1, 0]), Instruction(Ry(-ϕ/2), [0]), Instruction(CNot(), [1, 0]), Instruction(CNot(), [0, 1])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(SingleExcitationPlus(ϕ)))), [1, 0])]
        # state vector for SingleExcitationPlus (from PennyLane)
        state_vector = 0.5 * [exp(im*ϕ/2),  cos(ϕ/2) - sin(ϕ/2), cos(ϕ/2) + sin(ϕ/2), exp(im*ϕ/2)]
        probability_amplitudes =  0.25 * [1, (cos(ϕ/2) - sin(ϕ/2))^2, (cos(ϕ/2) + sin(ϕ/2))^2, 1]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(SingleExcitationPlus(ϕ)) == 2
        @test inv(SingleExcitationPlus(ϕ)) == SingleExcitationPlus(-ϕ)
        @test SingleExcitationPlus(ϕ) ^ 0 == SingleExcitationPlus(0.0)
        @test SingleExcitationPlus(ϕ) ^ 2 == SingleExcitationPlus(2*ϕ)
        @test SingleExcitationPlus(ϕ) ^ -1 == SingleExcitationPlus(-ϕ)
        @test SingleExcitationPlus(ϕ) ^ - 3 == SingleExcitationPlus(-3*ϕ)
    end
	
    @testset "Single excitation minus" begin
        ϕ  = 3.56
        nq = 2
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(SingleExcitationMinus(ϕ), [0, 1])]
        # instructions for the gate decomposition (from PennyLane)
        de_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(X(), [0]), Instruction(X(), [1]), Instruction(CPhaseShift(-ϕ/2), [1, 0]), Instruction(X(), [0]), Instruction(X(), [1]), Instruction(CPhaseShift(-ϕ/2), [0, 1]), Instruction(CNot(), [0, 1]), Instruction(Ry(ϕ/2), [0]), Instruction(CNot(), [1, 0]), Instruction(Ry(-ϕ/2), [0]), Instruction(CNot(), [1, 0]), Instruction(CNot(), [0, 1])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(SingleExcitationMinus(ϕ)))), [1, 0])]
        # state vector for SingleExcitationMinus (from PennyLane)
        state_vector = 0.5 * [exp(-im*ϕ/2),  cos(ϕ/2) - sin(ϕ/2), cos(ϕ/2) + sin(ϕ/2), exp(-im*ϕ/2)]
        probability_amplitudes =  0.25 * [1, (cos(ϕ/2) - sin(ϕ/2))^2, (cos(ϕ/2) + sin(ϕ/2))^2, 1]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("decomp", de_instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(SingleExcitationMinus(ϕ)) == 2
        @test inv(SingleExcitationMinus(ϕ)) == SingleExcitationMinus(-ϕ)
        @test SingleExcitationMinus(ϕ) ^ 0 == SingleExcitationMinus(0.0)
        @test SingleExcitationMinus(ϕ) ^ 2 == SingleExcitationMinus(2*ϕ)
        @test SingleExcitationMinus(ϕ) ^ -1 == SingleExcitationMinus(-ϕ)
        @test SingleExcitationMinus(ϕ) ^ - 3 == SingleExcitationMinus(-3*ϕ)
    end

    @testset "Double excitation minus" begin
        ϕ  = 3.56
        nq = 4
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(DoubleExcitationMinus(ϕ), [0, 1, 2, 3])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(DoubleExcitationMinus(ϕ)))), [0, 1, 2, 3])]
        state_vector =  0.5 * [exp(-im*ϕ/2.0), 0, 0, -sin(ϕ/2), exp(-im*ϕ/2.0), 0, 0, 0, exp(-im*ϕ/2.0), 0, 0, 0, cos(ϕ/2), 0, 0, 0]
        probability_amplitudes = 0.25 * [1, 0, 0, (-sin(ϕ/2))^2, 1, 0, 0, 0, 1, 0, 0, 0, (cos(ϕ/2))^2, 0, 0, 0]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(DoubleExcitationMinus(ϕ)) == 4
        @test inv(DoubleExcitationMinus(ϕ)) == DoubleExcitationMinus(-ϕ)
        @test DoubleExcitationMinus(ϕ) ^ 0 == DoubleExcitationMinus(0.0)
        @test DoubleExcitationMinus(ϕ) ^ 2 == DoubleExcitationMinus(2*ϕ)
        @test DoubleExcitationMinus(ϕ) ^ -1 == DoubleExcitationMinus(-ϕ)
        @test DoubleExcitationMinus(ϕ) ^ - 3 == DoubleExcitationMinus(-3*ϕ)
    end

    @testset "Double excitation plus" begin
        ϕ  = 3.56
        nq = 4
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(DoubleExcitationPlus(ϕ), [0, 1, 2, 3])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(DoubleExcitationPlus(ϕ)))), [0, 1, 2, 3])]
        state_vector =  0.5 * [exp(im*ϕ/2.0), 0, 0, -sin(ϕ/2), exp(im*ϕ/2.0), 0, 0, 0, exp(im*ϕ/2.0), 0, 0, 0, cos(ϕ/2), 0, 0, 0]
        probability_amplitudes = 0.25 * [1, 0, 0, (-sin(ϕ/2))^2, 1, 0, 0, 0, 1, 0, 0, 0, (cos(ϕ/2))^2, 0, 0, 0]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(DoubleExcitationPlus(ϕ)) == 4
        @test inv(DoubleExcitationPlus(ϕ)) == DoubleExcitationPlus(-ϕ)
        @test DoubleExcitationPlus(ϕ) ^ 0 == DoubleExcitationPlus(0.0)
        @test DoubleExcitationPlus(ϕ) ^ 2 == DoubleExcitationPlus(2*ϕ)
        @test DoubleExcitationPlus(ϕ) ^ -1 == DoubleExcitationPlus(-ϕ)
        @test DoubleExcitationPlus(ϕ) ^ - 3 == DoubleExcitationPlus(-3*ϕ)
    end

    @testset "FermionicSWAP" begin
        ϕ  = 3.56
        nq = 2
        instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(FermionicSWAP(ϕ), [0, 1])]
        # instructions for the matrix representation
        u_instructions = [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(Unitary(Matrix(matrix_rep(FermionicSWAP(ϕ)))), [1, 0])]
        # state vector for FermionicSWAP (from PennyLane)
        state_vector = 0.5 * [1, exp(im*ϕ/2.0)*cos(ϕ / 2.0) - im*exp(im*ϕ/2.0)*sin(ϕ/2.0), 
                              - im*exp(im*ϕ/2.0)*sin(ϕ/2.0) + exp(im*ϕ/2.0)*cos(ϕ/2.0), exp(im * ϕ)]
        probability_amplitudes = 0.25 * [1, (exp(im*ϕ/2.0)*cos(ϕ / 2.0) - im*exp(im*ϕ/2.0)*sin(ϕ/2.0))^2, 
                              (-im*exp(im*ϕ/2.0)*sin(ϕ/2.0) + exp(im*ϕ/2.0)*cos(ϕ/2.0))^2, 1]
        @testset "Simulator $sim, instruction set $ix_label" for sim in (StateVectorSimulator, DensityMatrixSimulator),
            (ix_label, ixs) in (("raw", instructions), ("unitary", u_instructions))
            simulation = sim(nq, 0)
            simulation = evolve!(simulation, ixs)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
        @test qubit_count(FermionicSWAP(ϕ)) == 2
        @test inv(FermionicSWAP(ϕ)) == FermionicSWAP(-ϕ)
        @test FermionicSWAP(ϕ) ^ 0 == FermionicSWAP(0.0)
        @test FermionicSWAP(ϕ) ^ 2 == FermionicSWAP(2*ϕ)
        @test FermionicSWAP(ϕ) ^ -1 == FermionicSWAP(-ϕ)
        @test FermionicSWAP(ϕ) ^ - 3 == FermionicSWAP(-3*ϕ)
    end
end
