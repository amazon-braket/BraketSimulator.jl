using Test, BraketSimulator

@testset "IR validation" begin
    sim = DensityMatrixSimulator(2, 0)
    results = [(type="NotARealResult",)]
    @test_throws ErrorException BraketSimulator._validate_ir_results_compatibility(sim, results, Val(:OpenQASM))
    @test_throws ErrorException BraketSimulator._validate_ir_results_compatibility(sim, results, Val(:JAQCD))
    c = BraketSimulator.Circuit()
    BraketSimulator.add_instruction!(c, BraketSimulator.Instruction(BraketSimulator.Rx(BraketSimulator.FreeParameter(:α)), 0))
    @test_throws ErrorException BraketSimulator._validate_input_provided(c)
    c = BraketSimulator.Circuit()
    BraketSimulator.add_instruction!(c, BraketSimulator.Instruction(BraketSimulator.Rx(0.1), 0))
    @test isnothing(BraketSimulator._validate_input_provided(c))

    
    sim = DensityMatrixSimulator(2, 0)
    @test_logs (:warn, "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience.") BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:OpenQASM))
    @test_logs (:warn, "You are running a noise-free circuit on the density matrix simulator. Consider running this circuit on the state vector simulator: LocalSimulator(\"braket_sv_v2\") for a better user experience.") BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:JAQCD))
    sim = StateVectorSimulator(2, 0)
    BraketSimulator.add_instruction!(c, BraketSimulator.Instruction(BraketSimulator.BitFlip(0.1), 0))
    @test_throws ValidationError BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:OpenQASM))
    @test_throws ValidationError BraketSimulator._validate_ir_instructions_compatibility(sim, c, Val(:JAQCD))
end
