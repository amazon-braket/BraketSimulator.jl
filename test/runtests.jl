using Test

@testset "BraketSimulator Tests" begin
    # Run the original branched simulator tests
    include("test_branched_simulator.jl")

    # Run the new end-to-end tests
    include("test_branched_simulator_end_to_end.jl")
    
    # Run the branched simulator operators with OpenQASM tests
    include("test_branched_simulator_operators_openqasm.jl")
end
