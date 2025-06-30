using BraketSimulator

# Test dynamic qubit allocation in the branched simulator
function test_dynamic_qubit_allocation()
    println("Testing dynamic qubit allocation...")
    
    # Create a simulator with 100 shots
    # The 2 initial qubits is irrelevant as the branched simulator assigns qubits dynamically
    simulator = StateVectorSimulator(2, 100)
    
    qasm_source = """
			OPENQASM 3.0;
			bit b;
			qubit[2] q;
			int[32] count = 0;

			// Initialize qubit 0 to |0⟩
			// Keep measuring and flipping until we get a 1
			b = 0;
			while (b == 0 && count < 3) {
				h q[0];       // Put qubit 0 in superposition
				b = measure q[0];  // Measure qubit 0
				count = count + 1;
			}

			// Apply X to qubit 1 if we got a 1 within 3 attempts
			if (b == 1) {
				x q[1];
			}
			"""
        
    # Evolve the program using the branched simulator
    branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
    println(branched_sim.instruction_sequences)
    println(branched_sim.measurements)
    println(branched_sim.variables)
    state_00 = BraketSimulator.calculate_current_state(branched_sim, 1)
    state_01 = BraketSimulator.calculate_current_state(branched_sim, 2)
    state_10 = BraketSimulator.calculate_current_state(branched_sim, 3)
    state_11 = BraketSimulator.calculate_current_state(branched_sim, 4)
    println(state_00)
    println(state_01)
    println(state_10)
    println(state_11)
    return branched_sim
end

# Run the test
branched_sim = test_dynamic_qubit_allocation()

# Print some information about the final state
println("\nFinal state information:")
println("Number of active paths: $(length(branched_sim.active_paths))")
println("Total number of qubits: $(branched_sim.n_qubits)")

# Print measurement results for the first active path
if !isempty(branched_sim.active_paths)
    path_idx = first(branched_sim.active_paths)
    println("\nMeasurement results for path $path_idx:")
    for (qubit, result) in sort(collect(branched_sim.measurements[path_idx]))
        println("Qubit $qubit: $result")
    end
end
