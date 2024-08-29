using Test, LinearAlgebra, Logging, BraketSimulator, DataStructures

@testset "Gate operations" begin
    nq = 4
    ϕ = 0.23
    θ = 0.14
    λ = -1.17
    q1_mat = [0 -im; im 0] 
    q2_mat = [1 0 0 0; 0 1 0 0 ; 0 0 0 -im; 0 0 im 0] 
    @testset "Inverted gates" begin
        @testset "1q Gate $g" for g in [BraketSimulator.X(),
                                        BraketSimulator.Y(),
                                        BraketSimulator.Z(),
                                        BraketSimulator.H(),
                                        BraketSimulator.GPi(ϕ),
                                        BraketSimulator.GPi2(ϕ),
                                        BraketSimulator.Rx(ϕ),
                                        BraketSimulator.Ry(ϕ),
                                        BraketSimulator.Rz(ϕ),
                                        BraketSimulator.PRx(θ, ϕ),
                                        BraketSimulator.PhaseShift(ϕ),
                                        BraketSimulator.S(),
                                        BraketSimulator.Si(),
                                        BraketSimulator.T(),
                                        BraketSimulator.Ti(),
                                        BraketSimulator.V(),
                                        BraketSimulator.Vi(),
                                        BraketSimulator.U(θ, ϕ, λ),
                                        BraketSimulator.Unitary(q1_mat),
                                        BraketSimulator.MultiQubitPhaseShift{1}(ϕ),
                                       ]
            @test inv(BraketSimulator.matrix_rep(g)) ≈ BraketSimulator.matrix_rep(inv(g))
            @test BraketSimulator.matrix_rep(inv(g)) * BraketSimulator.matrix_rep(g) ≈ Diagonal(ones(2))
            sim = StateVectorSimulator(nq, 0)
            instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1]
            sim = BraketSimulator.evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = BraketSimulator.evolve!(new_sim, [BraketSimulator.Instruction(g, [0]), BraketSimulator.Instruction(inv(g), [0])])
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "2q Gate $g" for g in [BraketSimulator.XX(ϕ),
                                        BraketSimulator.YY(ϕ),
                                        BraketSimulator.ZZ(ϕ),
                                        BraketSimulator.XY(ϕ),
                                        BraketSimulator.ECR(),
                                        BraketSimulator.Swap(),
                                        BraketSimulator.ISwap(),
                                        BraketSimulator.PSwap(ϕ),
                                        BraketSimulator.MS(θ, ϕ, λ),
                                        BraketSimulator.CPhaseShift(ϕ),
                                        BraketSimulator.CPhaseShift00(ϕ),
                                        BraketSimulator.CPhaseShift01(ϕ),
                                        BraketSimulator.CPhaseShift10(ϕ),
                                        BraketSimulator.SingleExcitation(ϕ),
                                        BraketSimulator.Unitary(q2_mat),
                                        BraketSimulator.CNot(),
                                        BraketSimulator.CY(),
                                        BraketSimulator.CZ(),
                                        BraketSimulator.CV(),
                                        BraketSimulator.MultiQubitPhaseShift{2}(ϕ),
                                        BraketSimulator.Control(BraketSimulator.MultiQubitPhaseShift{2}(ϕ), (0,0)),
                                        BraketSimulator.Control(BraketSimulator.X(), (1,)),
                                        BraketSimulator.Control(BraketSimulator.U(θ, ϕ, λ), (1,))]
            @test inv(BraketSimulator.matrix_rep(g)) ≈ BraketSimulator.matrix_rep(inv(g))
            @test BraketSimulator.matrix_rep(inv(g)) * BraketSimulator.matrix_rep(g) ≈ Diagonal(ones(4))
            sim = StateVectorSimulator(nq, 0)
            instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1]
            sim = BraketSimulator.evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = BraketSimulator.evolve!(new_sim, [BraketSimulator.Instruction(g, [0, 1]), BraketSimulator.Instruction(inv(g), [0, 1])])
            @test new_sim.state_vector ≈ sim.state_vector atol=1e-5
        end
        @testset "3q Gate $g" for g in [BraketSimulator.CCNot(), BraketSimulator.CSwap()]
            @test inv(BraketSimulator.matrix_rep(g)) ≈ BraketSimulator.matrix_rep(inv(g))
            @test BraketSimulator.matrix_rep(inv(g)) * BraketSimulator.matrix_rep(g) ≈ Diagonal(ones(8))
            sim = StateVectorSimulator(nq, 0)
            instructions = [BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1]
            sim = BraketSimulator.evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = BraketSimulator.evolve!(new_sim, [BraketSimulator.Instruction(g, [0, 1, 3]), BraketSimulator.Instruction(inv(g), [0, 1, 3])])
            @test new_sim.state_vector ≈ sim.state_vector
        end
    end
    @testset "Exponentiated gates - pow $pow" for pow in (-0.5, 0, 0.5, 1, 2, 2.5, 3, 4, -2)
        @testset "1q Gate $g" for g in [BraketSimulator.X(),
                                        BraketSimulator.Y(),
                                        BraketSimulator.Z(),
                                        BraketSimulator.H(),
                                        BraketSimulator.GPi(ϕ),
                                        BraketSimulator.GPi2(ϕ),
                                        BraketSimulator.Rx(ϕ),
                                        BraketSimulator.Ry(ϕ),
                                        BraketSimulator.Rz(ϕ),
                                        BraketSimulator.PRx(θ, ϕ),
                                        BraketSimulator.PhaseShift(ϕ),
                                        BraketSimulator.S(),
                                        BraketSimulator.Si(),
                                        BraketSimulator.T(),
                                        BraketSimulator.Ti(),
                                        BraketSimulator.V(),
                                        BraketSimulator.Vi(),
                                        BraketSimulator.U(θ, ϕ, λ),
                                        BraketSimulator.Unitary(q1_mat),
                                        BraketSimulator.MultiQubitPhaseShift{1}(ϕ)
                                       ]
            sim = StateVectorSimulator(nq, 0)
            new_sim = StateVectorSimulator(nq, 0)
            instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(g^pow, [0])])
            new_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(g)^pow)), [0])])
            @test BraketSimulator.matrix_rep(g)^pow ≈ BraketSimulator.matrix_rep(g^pow)
            sim = BraketSimulator.evolve!(sim, instructions)
            new_sim = BraketSimulator.evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "2q Gate $g" for g in [BraketSimulator.XX(ϕ),
                                        BraketSimulator.YY(ϕ),
                                        BraketSimulator.ZZ(ϕ),
                                        BraketSimulator.XY(ϕ),
                                        BraketSimulator.ECR(),
                                        BraketSimulator.Swap(),
                                        BraketSimulator.ISwap(),
                                        BraketSimulator.PSwap(ϕ),
                                        BraketSimulator.MS(θ, ϕ, λ),
                                        BraketSimulator.CPhaseShift(ϕ),
                                        BraketSimulator.CPhaseShift00(ϕ),
                                        BraketSimulator.CPhaseShift01(ϕ),
                                        BraketSimulator.CPhaseShift10(ϕ),
                                        BraketSimulator.SingleExcitation(ϕ),
                                        BraketSimulator.Unitary(q2_mat),
                                        BraketSimulator.Control(BraketSimulator.MultiQubitPhaseShift{2}(ϕ), (0,0)),
                                        BraketSimulator.MultiQubitPhaseShift{2}(ϕ),]
            sim              = StateVectorSimulator(nq, 0)
            new_sim          = StateVectorSimulator(nq, 0)
            new_ixs          = [BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(g)^pow)), [0, 1])]
            instructions     = vcat([BraketSimulator.Instruction(BraketSimulator.H(), q) for q in 0:nq-1], [BraketSimulator.Instruction(g^pow, [0, 1])])
            new_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), q) for q in 0:nq-1], new_ixs)
            @test Matrix(BraketSimulator.matrix_rep(g))^pow ≈ Matrix(BraketSimulator.matrix_rep(g^pow))
            sim     = BraketSimulator.evolve!(sim, instructions)
            new_sim = BraketSimulator.evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "2q Gate $g" for (g,tg) in [(BraketSimulator.CNot(), BraketSimulator.X()),
                                             (BraketSimulator.CY(), BraketSimulator.Y()),
                                             (BraketSimulator.CZ(), BraketSimulator.Z()),
                                             (BraketSimulator.CV(), BraketSimulator.V()),
                                             (BraketSimulator.Control(BraketSimulator.X(), (1,)), BraketSimulator.X()),
                                             (BraketSimulator.Control(BraketSimulator.U(θ, ϕ, λ), (1,)), BraketSimulator.U(θ, ϕ, λ))]
            sim              = StateVectorSimulator(nq, 0)
            new_sim          = StateVectorSimulator(nq, 0)
            new_ixs          = [BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(tg)^pow)), (1,)), [0, 1])]
            instructions     = vcat([BraketSimulator.Instruction(BraketSimulator.H(), q) for q in 0:nq-1], [BraketSimulator.Instruction(g^pow, [0, 1])])
            new_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), q) for q in 0:nq-1], new_ixs)
            sim     = BraketSimulator.evolve!(sim, instructions)
            new_sim = BraketSimulator.evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "3q Gate $g" for g in [BraketSimulator.CCNot(), BraketSimulator.CSwap()]
            nq  = 5
            sim = StateVectorSimulator(nq, 0)
            new_sim = StateVectorSimulator(nq, 0)
            instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(g^pow, [0, 1, 3])])
            new_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(g)^pow)), [0, 1, 3])])
            if !isa(g^pow, BraketSimulator.I)
                @test BraketSimulator.matrix_rep(g)^pow ≈ BraketSimulator.matrix_rep(g^pow)
            end
            sim = BraketSimulator.evolve!(sim, instructions)
            new_sim = BraketSimulator.evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "4q Gate $g" for g in [BraketSimulator.DoubleExcitation(ϕ),]
            nq  = 5
            sim = StateVectorSimulator(nq, 0)
            new_sim = StateVectorSimulator(nq, 0)
            instructions     = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(g^pow, [0, 1, 3, 2])])
            new_instructions = vcat([BraketSimulator.Instruction(BraketSimulator.H(), [q]) for q in 0:nq-1], [BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(BraketSimulator.matrix_rep(g)))^pow, [0, 1, 3, 2])])
            if !isa(g^pow, BraketSimulator.I)
                @test BraketSimulator.matrix_rep(g)^pow ≈ BraketSimulator.matrix_rep(g^pow)
            end
            sim     = BraketSimulator.evolve!(sim, instructions)
            new_sim = BraketSimulator.evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
    end
end
