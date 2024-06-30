using Test, Braket, BraketSimulator

@testset "IR validation" begin
    sim = DensityMatrixSimulator(2, 0)
    results = [(type="NotARealResult",)]
    @test_throws ErrorException BraketSimulator._validate_ir_results_compatibility(sim, results, Val(:OpenQASM))
    @test_throws ErrorException BraketSimulator._validate_ir_results_compatibility(sim, results, Val(:JAQCD))
    
    @test_throws ErrorException BraketSimulator._validate_input_provided(Circuit([(Rx(FreeParameter(:α)), 0)]))
    @test isnothing(BraketSimulator._validate_input_provided(Circuit([(Rx(0.1), 0)])))

    
    sim = DensityMatrixSimulator(2, 0)
    c = Circuit([(Rx(0.1), 0)])
    @test_logs (:warn, "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience.") BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:OpenQASM))
    @test_logs (:warn, "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience.") BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:JAQCD))
    sim = StateVectorSimulator(2, 0)
    c = Circuit([(BitFlip(0.1), 0)])
    @test_throws ValidationError BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:OpenQASM))
    @test_throws ValidationError BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:JAQCD))
end
