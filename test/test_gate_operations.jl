using Test, LinearAlgebra, Logging, Braket, BraketSimulator, DataStructures

using Braket: Instruction
using BraketSimulator: DoubleExcitation, SingleExcitation, matrix_rep, Control

@testset "Gate operations" begin
    nq = 4
    ϕ = 0.23
    θ = 0.14
    λ = -1.17
    q1_mat = [0 -im; im 0] 
    q2_mat = [1 0 0 0; 0 1 0 0 ; 0 0 0 -im; 0 0 im 0] 
    @testset "Inverted gates" begin
        @testset for g in [X(), Y(), Z(), H(), GPi(ϕ), GPi2(ϕ), Rx(ϕ), Ry(ϕ), Rz(ϕ), PhaseShift(ϕ), S(), Si(), T(), Ti(), V(), Vi(), U(θ, ϕ, λ), Unitary(q1_mat)]
            @test inv(matrix_rep(g)) ≈ matrix_rep(inv(g))
            @test matrix_rep(inv(g)) * matrix_rep(g) ≈ Diagonal(ones(2))
            sim = StateVectorSimulator(nq, 0)
            instructions = [Instruction(H(), [q]) for q in 0:nq-1]
            sim = evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = evolve!(new_sim, [Instruction(g, [0]), Instruction(inv(g), [0])])
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset for g in [XX(ϕ), YY(ϕ), ZZ(ϕ), XY(ϕ), ECR(), Swap(), ISwap(), PSwap(ϕ), MS(θ, ϕ, λ), CPhaseShift(ϕ), CPhaseShift00(ϕ), CPhaseShift01(ϕ), CPhaseShift10(ϕ), Unitary(q2_mat), CNot(), CY(), CZ(), CV()]
            @test inv(matrix_rep(g)) ≈ matrix_rep(inv(g))
            @test matrix_rep(inv(g)) * matrix_rep(g) ≈ Diagonal(ones(4))
            sim = StateVectorSimulator(nq, 0)
            instructions = [Instruction(H(), [q]) for q in 0:nq-1]
            sim = evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = evolve!(new_sim, [Instruction(g, [0, 1]), Instruction(inv(g), [0, 1])])
            @test new_sim.state_vector ≈ sim.state_vector atol=1e-5
        end
        @testset for g in [CCNot(), CSwap()]
            @test inv(matrix_rep(g)) ≈ matrix_rep(inv(g))
            @test matrix_rep(inv(g)) * matrix_rep(g) ≈ Diagonal(ones(8))
            sim = StateVectorSimulator(nq, 0)
            instructions = [Instruction(H(), [q]) for q in 0:nq-1]
            sim = evolve!(sim, instructions)
            new_sim = copy(sim)
            new_sim = evolve!(new_sim, [Instruction(g, [0, 1, 3]), Instruction(inv(g), [0, 1, 3])])
            @test new_sim.state_vector ≈ sim.state_vector
        end
    end
    @testset "Exponentiated gates" begin
        @testset "Gate $g, pow $pow" for g in [X(), Y(), Z(), H(), GPi(ϕ), GPi2(ϕ), Rx(ϕ), Ry(ϕ), Rz(ϕ), PhaseShift(ϕ), S(), Si(), T(), Ti(), V(), Vi(), U(θ, ϕ, λ), Unitary(q1_mat), MultiQubitPhaseShift{1}(ϕ)], pow in (0, 1, 2, 3, 4)
            sim = StateVectorSimulator(nq, 0)
            new_sim = StateVectorSimulator(nq, 0)
            instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(g^pow, [0])])
            new_instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(Unitary(Matrix(matrix_rep(g)^pow)), [0])])
            @test matrix_rep(g)^pow ≈ matrix_rep(g^pow)
            sim = evolve!(sim, instructions)
            new_sim = evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "Gate $g, pow $pow" for g in [XX(ϕ), YY(ϕ), ZZ(ϕ), XY(ϕ), CNot(), CY(), CZ(), CV(), ECR(), Swap(), ISwap(), PSwap(ϕ), MS(θ, ϕ, λ), CPhaseShift(ϕ), CPhaseShift00(ϕ), CPhaseShift01(ϕ), CPhaseShift10(ϕ), Unitary(q2_mat), MultiQubitPhaseShift{2}(ϕ), Control(MultiQubitPhaseShift{2}(ϕ), (0,0)), Control(X(), (1,))], pow in (0, 1, 2, 3, 4)
            sim              = StateVectorSimulator(nq, 0)
            new_sim          = StateVectorSimulator(nq, 0)
            instructions     = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(g^pow, [0, 1])])
            new_ixs          = [Instruction(g, [0, 1]) for p in 1:pow]
            new_instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], new_ixs)
            if !isa(g^pow, Braket.I)
                @test Matrix(matrix_rep(g))^pow ≈ Matrix(matrix_rep(g^pow))
            end
            sim     = evolve!(sim, instructions)
            new_sim = evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
        @testset "Gate $g, pow $pow" for g in [CCNot(), CSwap()], pow in (0, 1, 2, 3, 4)
            sim = StateVectorSimulator(nq, 0)
            new_sim = StateVectorSimulator(nq, 0)
            instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(g^pow, [0, 1, 3])])
            new_instructions = vcat([Instruction(H(), [q]) for q in 0:nq-1], [Instruction(Unitary(Matrix(matrix_rep(g)^pow)), [0, 1, 3])])
            if !isa(g^pow, Braket.I)
                @test matrix_rep(g)^pow ≈ matrix_rep(g^pow)
            end
            sim = evolve!(sim, instructions)
            new_sim = evolve!(new_sim, new_instructions)
            @test new_sim.state_vector ≈ sim.state_vector
        end
    end
end
