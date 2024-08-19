using BraketSimulator, Test

@testset "QubitSet and Qubit" begin
    @testset "ctors" begin
        @test BraketSimulator.Qubit(1) == BraketSimulator.Qubit(1.0)
        @test BraketSimulator.Qubit(BigFloat(1.0)) == BraketSimulator.Qubit(1.0)
        @test BraketSimulator.Qubit(BraketSimulator.Qubit(1)) == BraketSimulator.Qubit(1.0)
        @test BraketSimulator.QubitSet(nothing) == BraketSimulator.QubitSet()
    end
    @testset "equality" begin
        @test BraketSimulator.Qubit(1) == Int8(1)
        @test BigInt(10) == BraketSimulator.Qubit(10)
        @test 10 == BraketSimulator.Qubit(10)
        @test BraketSimulator.Qubit(10) == BigInt(10)
    end
    @testset "Convert to Int" begin
        @test Int(BraketSimulator.Qubit(1)) == 1
        @test convert(Int, BraketSimulator.Qubit(1)) == 1
    end
    @testset "show" begin
        s = sprint(show, BraketSimulator.Qubit(1))
        @test s == "Qubit(1)"
    end
    @testset "QubitSet indexing" begin
        qs = BraketSimulator.QubitSet(0, 1, 2)
        @test copy(qs) == qs 
        @test length(qs) == 3
        @test lastindex(qs) == 3
        @test 3 ∉ qs
        @test !isempty(qs)
        o = popfirst!(qs)
        @test o == 0
        @test length(qs) == 2
        q1 = BraketSimulator.QubitSet(0, 1, 2)
        q2 = BraketSimulator.QubitSet(2, 4, 3)
        @test q1[1:2] == BraketSimulator.QubitSet(0, 1)
        @test intersect(q1, q2) == BraketSimulator.QubitSet(2)
        @test sort(q2) == BraketSimulator.QubitSet(2, 3, 4)
    end
end
