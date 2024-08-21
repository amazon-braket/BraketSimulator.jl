using Test, BraketSimulator, LinearAlgebra

CCNot_mat = round.(reduce(hcat, [[1.0, 0, 0, 0, 0, 0, 0, 0],
             [0, 1.0, 0, 0, 0, 0, 0, 0],
             [0, 0, 1.0, 0, 0, 0, 0, 0],
             [0, 0, 0, 1.0, 0, 0, 0, 0],
             [0, 0, 0, 0, 1.0, 0, 0, 0],
             [0, 0, 0, 0, 0, 1.0, 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 1.0],
             [0, 0, 0, 0, 0, 0, 1.0, 0]]), digits=8)
ECR_mat = round.(reduce(hcat, [[0, 0, 0.70710678, 0.70710678im],
                        [0, 0, 0.70710678im, 0.70710678],
                        [0.70710678, -0.70710678im, 0, 0],
                        [-0.70710678im, 0.70710678, 0, 0]]), digits=8)
T_mat = round.(reduce(hcat, [[1.0, 0], [0, 0.70710678 + 0.70710678im]]), digits=8)

@testset "Gates" begin
    @testset for g in (BraketSimulator.H(),
                       BraketSimulator.I(),
                       BraketSimulator.X(),
                       BraketSimulator.Y(),
                       BraketSimulator.Z(),
                       BraketSimulator.S(),
                       BraketSimulator.Si(),
                       BraketSimulator.T(),
                       BraketSimulator.Ti(),
                       BraketSimulator.V(),
                       BraketSimulator.Vi())
        @test qubit_count(g) == 1
        @test qubit_count(typeof(g)) == 1
        c = BraketSimulator.Circuit()
        BraketSimulator.add_instruction!(c, BraketSimulator.Instruction(g, 0))
        @test c.instructions == [BraketSimulator.Instruction(g, 0)]
        @test BraketSimulator.Parametrizable(g) == BraketSimulator.NonParametrized()
        @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
        @test BraketSimulator.n_angles(g) == 0
        @test BraketSimulator.angles(g) == ()
        @test BraketSimulator.qubit_count(g) == 1
    end
    @testset for g in (BraketSimulator.CNot(),
                       BraketSimulator.Swap(),
                       BraketSimulator.ISwap(),
                       BraketSimulator.CV(),
                       BraketSimulator.CY(),
                       BraketSimulator.CZ(),
                       BraketSimulator.ECR())
        @test qubit_count(g) == 2
        @test qubit_count(typeof(g)) == 2
        ix = BraketSimulator.Instruction(g, [0, 1])
        c = BraketSimulator.Circuit()
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        ix = BraketSimulator.Instruction(g, [10, 1])
        c = BraketSimulator.Circuit()
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        @test BraketSimulator.Parametrizable(g) == BraketSimulator.NonParametrized()
        @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
        @test BraketSimulator.n_angles(g) == 0
        @test BraketSimulator.qubit_count(g) == 2
    end
    @testset for g in (BraketSimulator.CCNot(), BraketSimulator.CSwap())
        @test qubit_count(g) == 3
        @test qubit_count(typeof(g)) == 3
        c = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(g, [0, 1, 2])
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        @test BraketSimulator.Parametrizable(g) == BraketSimulator.NonParametrized()
        @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
        @test BraketSimulator.n_angles(g) == 0
        @test BraketSimulator.qubit_count(g) == 3
    end
    α = BraketSimulator.FreeParameter(:α)
    β = BraketSimulator.FreeParameter(:β)
    γ = BraketSimulator.FreeParameter(:γ)
    @test copy(α) === α
    @testset for (angle1, angle2, angle3) in ((rand(), rand(), rand()), (π, rand(), rand()))
        @testset for g in (BraketSimulator.Rx(angle1), BraketSimulator.Ry(angle1), BraketSimulator.Rz(angle1), BraketSimulator.PhaseShift(angle1))
            @test qubit_count(g) == 1
            @test qubit_count(typeof(g)) == 1
            ix = BraketSimulator.Instruction(g, 0)
            c = BraketSimulator.Circuit()
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
            @test BraketSimulator.n_angles(g) == 1
            @test BraketSimulator.qubit_count(g) == 1
        end
        @testset for g in (BraketSimulator.Rx, BraketSimulator.Ry, BraketSimulator.Rz, BraketSimulator.PhaseShift)
            @test qubit_count(g) == 1
            ix = BraketSimulator.Instruction(g(α), 0)
            c = BraketSimulator.Circuit()
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g(α)) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g(α)) == BraketSimulator.FreeParameter[α]
            bound_g = BraketSimulator.bind_value!(BraketSimulator.Parametrized(), g(α), Dict(:α=>angle1))
            @test bound_g == g(angle1)
        end
        @testset for g in (BraketSimulator.PSwap(angle1),
                           BraketSimulator.XY(angle1),
                           BraketSimulator.CPhaseShift(angle1),
                           BraketSimulator.CPhaseShift00(angle1),
                           BraketSimulator.CPhaseShift01(angle1),
                           BraketSimulator.CPhaseShift10(angle1),
                           BraketSimulator.XX(angle1),
                           BraketSimulator.YY(angle1),
                           BraketSimulator.ZZ(angle1))
            @test qubit_count(g) == 2
            @test qubit_count(typeof(g)) == 2
            ix = BraketSimulator.Instruction(g, [0, 1])
            c = BraketSimulator.Circuit()
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
            @test BraketSimulator.n_angles(g) == 1
            @test BraketSimulator.qubit_count(g) == 2
        end
        @testset for g in (BraketSimulator.PSwap, BraketSimulator.XY, BraketSimulator.CPhaseShift, BraketSimulator.CPhaseShift00, BraketSimulator.CPhaseShift01, BraketSimulator.CPhaseShift10, BraketSimulator.XX, BraketSimulator.YY, BraketSimulator.ZZ)
            @test qubit_count(g) == 2
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α), [0, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α), [10, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g(α)) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g(α)) == BraketSimulator.FreeParameter[α]
            bound_g = BraketSimulator.bind_value!(BraketSimulator.Parametrized(), g(α), Dict(:α=>angle1))
            @test bound_g == g(angle1)
        end
        @testset for g in (BraketSimulator.MS(angle1, angle2, angle3),)
            @test qubit_count(g) == 2
            c  = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g, [0, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            c  = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g, [10, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
            @test BraketSimulator.n_angles(g) == 3
            @test BraketSimulator.qubit_count(g) == 2
        end
        @testset for g in (BraketSimulator.MS,)
            @test qubit_count(g) == 2
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α, β, γ), [0, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α, β, γ), [10, 1])
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g(α, β, γ)) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g(α, β, γ)) == BraketSimulator.FreeParameter[α, β, γ]
            bound_g = BraketSimulator.bind_value!(BraketSimulator.Parametrized(), g(α, β, γ), Dict(:α=>angle1, :β=>angle2, :γ=>angle3))
            @test bound_g == g(angle1, angle2, angle3)
        end
        @testset for g in (BraketSimulator.U(angle1, angle2, angle3),)
            @test qubit_count(g) == 1
            c  = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g, 0)
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
            @test BraketSimulator.n_angles(g) == 3
            @test BraketSimulator.qubit_count(g) == 1
        end
        @testset for g in (BraketSimulator.U,)
            @test qubit_count(g) == 1
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α, β, γ), 0)
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g(α, β, γ)) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g(α, β, γ)) == BraketSimulator.FreeParameter[α, β, γ]
            bound_g = BraketSimulator.bind_value!(BraketSimulator.Parametrized(), g(α, β, γ), Dict(:α=>angle1, :β=>angle2, :γ=>angle3))
            @test bound_g == g(angle1, angle2, angle3)
        end
        @testset for g in (BraketSimulator.PRx(angle1, angle2),)
            c  = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g, 0)
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g) == BraketSimulator.FreeParameter[]
            @test BraketSimulator.n_angles(g) == 2
            @test BraketSimulator.qubit_count(g) == 1
        end
        @testset for g in (BraketSimulator.PRx,)
            @test qubit_count(g) == 1
            c = BraketSimulator.Circuit()
            ix = BraketSimulator.Instruction(g(α, β), 0)
            BraketSimulator.add_instruction!(c, ix)
            @test c.instructions == [ix]
            @test BraketSimulator.Parametrizable(g(α, β)) == BraketSimulator.Parametrized()
            @test BraketSimulator.parameters(g(α, β)) == BraketSimulator.FreeParameter[α, β]
            bound_g = BraketSimulator.bind_value!(BraketSimulator.Parametrized(), g(α, β), Dict(:α=>angle1, :β=>angle2))
            @test bound_g == g(angle1, angle2)
        end
    end
    @testset "g = Unitary" begin
        X = complex([0. 1.; 1. 0.])
        Y = [0. -im; im 0.]
        g = BraketSimulator.Unitary(X)
        c = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(g, 0)
        @test BraketSimulator.qubit_count(g) == 1
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        g = BraketSimulator.Unitary(kron(X, X))
        @test BraketSimulator.qubit_count(g) == 2
        c = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(g, [0, 1])
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        c  = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(g, [0, 1])
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        g = BraketSimulator.Unitary(kron(Y, Y))
        @test BraketSimulator.qubit_count(g) == 2
        ix = BraketSimulator.Instruction(g, [0, 1])
        c = BraketSimulator.Circuit()
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
    end
    @testset "Angled 3 qubit gates" begin
        # build some "fake" (for now) 3 qubit gates to test gate applicators
        mutable struct CCPhaseShift <: BraketSimulator.AngledGate{1}
            angle::NTuple{1, Union{Real, BraketSimulator.FreeParameter}}
            pow_exponent::Float64
            CCPhaseShift(angle::T, pow_exponent::Float64=1.0) where {T<:NTuple{1, Union{Real, FreeParameter}}} = new(angle, pow_exponent)
        end
        mutable struct CXX <: BraketSimulator.AngledGate{1}
            angle::NTuple{1, Union{Float64, BraketSimulator.FreeParameter}}
            pow_exponent::Float64
            CXX(angle::T, pow_exponent::Float64=1.0) where {T<:NTuple{1, Union{Real, FreeParameter}}} = new(angle, pow_exponent)
        end
        angle = rand()
        c  = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(CCPhaseShift(angle), [0, 1, 2])
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
        c  = BraketSimulator.Circuit()
        ix = BraketSimulator.Instruction(CXX(angle), [0, 1, 2])
        BraketSimulator.add_instruction!(c, ix)
        @test c.instructions == [ix]
    end
end
