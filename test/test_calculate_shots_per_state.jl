using Test
using BraketSimulator

@testset "calculate_shots_per_state" begin
    # Create a simple circuit with a Hadamard gate and measurement
    qasm_code = """
    OPENQASM 3.0;
    bit b;
    qubit q;
    h q;
    b = measure q;
    """
    
    # Create the program
    program = OpenQasmProgram(
        braketSchemaHeader("braket.ir.openqasm.program", "1"),
        qasm_code,
        nothing
    )
    
    # Create a simulator with 1000 shots
    shots = 1000
    base_simulator = StateVectorSimulator(0, shots)
    simulator = BranchedSimulatorOperators(base_simulator)
    
    # Evolve the circuit
    branched_sim = evolve_branched_operators(simulator, new_to_circuit(qasm_code), Dict{String, Any}())
    
    # Calculate shots per state
    shots_per_state = calculate_shots_per_state(branched_sim)
    
    # We expect approximately 500 shots for state "0" and 500 for state "1"
    @test haskey(shots_per_state, "0")
    @test haskey(shots_per_state, "1")
    @test abs(shots_per_state["0"] - 500) < 50  # Allow some variance due to randomness
    @test abs(shots_per_state["1"] - 500) < 50
    @test shots_per_state["0"] + shots_per_state["1"] == shots  # Total should equal original shots
    
    # Test with a more complex circuit (2 qubits in superposition)
    qasm_code2 = """
    OPENQASM 3.0;
    bit[2] b;
    qubit[2] q;
    h q[0];
    h q[1];
    b[0] = measure q[0];
    b[1] = measure q[1];
    """
    
    # Create the program
    program2 = OpenQasmProgram(
        braketSchemaHeader("braket.ir.openqasm.program", "1"),
        qasm_code2,
        nothing
    )
    
    # Create a simulator with 1000 shots
    shots = 1000
    base_simulator = StateVectorSimulator(0, shots)
    simulator = BranchedSimulatorOperators(base_simulator)
    
    # Evolve the circuit
    branched_sim = evolve_branched_operators(simulator, new_to_circuit(qasm_code2), Dict{String, Any}())
    
    # Calculate shots per state
    shots_per_state = calculate_shots_per_state(branched_sim)
    
    # We expect approximately 250 shots for each of the 4 possible states
    @test haskey(shots_per_state, "00")
    @test haskey(shots_per_state, "01")
    @test haskey(shots_per_state, "10")
    @test haskey(shots_per_state, "11")
    
    # Each state should have approximately 250 shots
    for state in ["00", "01", "10", "11"]
        @test abs(shots_per_state[state] - 250) < 50
    end
    
    # Total should equal original shots
    @test sum(values(shots_per_state)) == shots
    
    # Print the results
    println("Shots per state for 2-qubit circuit:")
    for (state, count) in shots_per_state
        println("State $state: $count shots")
    end
end
