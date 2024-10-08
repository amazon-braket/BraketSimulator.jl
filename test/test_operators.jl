using BraketSimulator, BraketSimulator.Dates, Test

@testset "Operators" begin
    for op in (BraketSimulator.Measure(),
               BraketSimulator.Reset(),
               BraketSimulator.Barrier(),
               BraketSimulator.Delay(Microsecond(200))
              )
        @test BraketSimulator.qubit_count(op)         == 1
        @test BraketSimulator.qubit_count(typeof(op)) == 1
        @test BraketSimulator.parameters(op)          == BraketSimulator.FreeParameter[]
    end
end
