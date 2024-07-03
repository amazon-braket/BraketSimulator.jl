using Test, BraketSimulator, LinearAlgebra

struct CustomNoise <: BraketSimulator.Noise end

@testset "Noises" begin
    @testset for noise in (BraketSimulator.BitFlip, BraketSimulator.PhaseFlip, BraketSimulator.AmplitudeDamping, BraketSimulator.PhaseDamping, BraketSimulator.Depolarizing)
        n = noise(0.1)
        @test BraketSimulator.qubit_count(n) == 1
        ix = BraketSimulator.Instruction(n, 0)
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset "noise = PauliChannel" begin
        n = BraketSimulator.PauliChannel(0.1, 0.2, 0.1)
        @test BraketSimulator.qubit_count(n) == 1
        ix = BraketSimulator.Instruction(n, 0)
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset "noise = MultiQubitPauliChannel{1}" begin
        n = BraketSimulator.MultiQubitPauliChannel{1}(Dict("X"=>0.1, "Y"=>0.2))
        @test BraketSimulator.qubit_count(n) == 1
        ix = BraketSimulator.Instruction(n, 0)
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset "noise = BraketSimulator.TwoQubitPauliChannel" begin
        n = BraketSimulator.TwoQubitPauliChannel(Dict("XX"=>0.1, "YY"=>0.2))
        @test BraketSimulator.qubit_count(n) == 2
        ix = BraketSimulator.Instruction(n, [0, 1])
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset for noise in (BraketSimulator.TwoQubitDephasing, BraketSimulator.TwoQubitDepolarizing)
        n = noise(0.4)
        @test BraketSimulator.qubit_count(n) == 2
        ix = BraketSimulator.Instruction(n, [0, 1])
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset "noise = GeneralizedAmplitudeDamping" begin
        n = BraketSimulator.GeneralizedAmplitudeDamping(0.1, 0.2)
        @test BraketSimulator.qubit_count(n) == 1
        ix = BraketSimulator.Instruction(n, 0)
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.Parametrized()
    end
    @testset "noise = Kraus" begin
        mat = complex([0 1; 1 0])
        n = BraketSimulator.Kraus([mat])
        @test BraketSimulator.qubit_count(n) == 1
        ix = BraketSimulator.Instruction(n, 0)
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.NonParametrized()
    end
    @testset "fallback" begin
        n = CustomNoise()
        @test BraketSimulator.bind_value!(n, Dict(:theta=>0.1)) === n
        @test BraketSimulator.parameters(n) == BraketSimulator.FreeParameter[]
        @test BraketSimulator.Parametrizable(n) == BraketSimulator.NonParametrized()
    end
    @test BraketSimulator.StructTypes.StructType(Noise) == BraketSimulator.StructTypes.AbstractType()
end
