using BraketSimulator

# Test dynamic qubit allocation in the branched simulator
function test_dynamic_qubit_allocation()
    println("Testing dynamic qubit allocation...")
    
    # Create a simulator with 100 shots
    # The 2 initial qubits is irrelevant as the branched simulator assigns qubits dynamically
    simulator = StateVectorSimulator(2, 100)
    
    # Define a program that creates registers dynamically
    qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[3] q;

        h q[0];       // Put qubit 0 in superposition
        h q[1];       // Put qubit 1 in superposition
        b[0] = measure q[0];  // Measure qubit 0
        b[1] = measure q[1];  // Measure qubit 1


        if (b[0] == 1) {
            if (b[1] == 1){  // Both measured 1
                x q[2];    // Apply X to qubit 2
            } else{
                h q[2];    // Apply H to qubit 2
            }
        } else {
            if (b[1] == 1) {  // Only second qubit measured 1
                z q[2];    // Apply Z to qubit 2
            }
        }
        // If both measured 0, do nothing to qubit 2
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
