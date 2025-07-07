using Test, LinearAlgebra, Random
using BraketSimulator

@testset "Qubit Aliasing with let keyword" begin
    @testset "Basic qubit aliasing" begin
        # Create a simple OpenQASM program with qubit aliasing
        qasm_source = """
        OPENQASM 3.0;
        qubit[5] q;
        
        // Create an alias for a single qubit
        let single_qubit = q[1];
        
        // Create an alias for a range of qubits
        let myreg = q[1:4];
        
        // Use the aliases
        h single_qubit;
        x myreg[2];  // This should apply X to q[3]
        """

        # Create a simulator
        simulator = StateVectorSimulator(5, 1000)

        # Evolve the program using the branched simulator operators
        branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

        # Calculate the current state
        state = BraketSimulator.calculate_current_state(branched_sim, 1)
        
        # Debug: Print the state vector
        @info "State vector" state
        
        # Verify that H was applied to q[1]
        # The state should have equal probability for q[1] being 0 or 1
        # Extract the probabilities for q[1] = 0 and q[1] = 1
        prob_q1_0 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 3) == 0)
        prob_q1_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 3) != 0)
        
        @test isapprox(prob_q1_0, 0.5, atol=1e-10)
        @test isapprox(prob_q1_1, 0.5, atol=1e-10)
        
        # Verify that X was applied to q[3]
        # The state should have q[3] = 1
        # Extract the probability for q[3] = 1
        prob_q3_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 1) != 0)
        
        @test isapprox(prob_q3_1, 1.0, atol=1e-10)
    end

    @testset "Aliasing with concatenation" begin
        # Create an OpenQASM program with concatenated qubit aliases
        qasm_source = """
        OPENQASM 3.0;
        qubit[3] q1;
        qubit[2] q2;
        
        // Create aliases with concatenation
        let combined = q1[0:1] ++ q2;
        
        // Use the alias
        h combined[0];  // Should apply H to q1[0]
        x combined[2];  // Should apply X to q2[0]
        """

        # Create a simulator
        simulator = StateVectorSimulator(5, 1000)

        # Evolve the program using the branched simulator operators
        branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

        # Calculate the current state
        state = BraketSimulator.calculate_current_state(branched_sim, 1)
        
        # Verify that H was applied to q1[0]
        # Extract the probabilities for q1[0] = 0 and q1[0] = 1
        prob_q1_0_0 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 4) == 0)
        prob_q1_0_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 4) != 0)
        
        @test isapprox(prob_q1_0_0, 0.5, atol=1e-10)
        @test isapprox(prob_q1_0_1, 0.5, atol=1e-10)
        
        # Verify that X was applied to q2[0] (which is q1[0:1] ++ q2[0], so index 3)
        # Extract the probability for q2[0] = 1
        prob_q2_0_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 1) != 0)
        
        @test isapprox(prob_q2_0_1, 1.0, atol=1e-10)
    end

    @testset "Nested aliasing" begin
        # Create an OpenQASM program with nested qubit aliases
        qasm_source = """
        OPENQASM 3.0;
        qubit[5] q;
        
        // Create first level alias
        let reg1 = q[0:2];
        
        // Create nested alias
        let reg2 = reg1[1:2];
        
        // Use the nested alias
        h reg2[0];  // Should apply H to q[1]
        """

        # Create a simulator
        simulator = StateVectorSimulator(5, 1000)

        # Evolve the program using the branched simulator operators
        branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

        # Calculate the current state
        state = BraketSimulator.calculate_current_state(branched_sim, 1)
        
        # Verify that H was applied to q[1]
        # Extract the probabilities for q[1] = 0 and q[1] = 1
        prob_q1_0 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 3) == 0)
        prob_q1_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 3) != 0)
        
        @test isapprox(prob_q1_0, 0.5, atol=1e-10)
        @test isapprox(prob_q1_1, 0.5, atol=1e-10)
    end

    @testset "Aliasing in local scopes" begin
        # Create an OpenQASM program with aliases in local scopes
        qasm_source = """
        OPENQASM 3.0;
        qubit[5] q;
        bit b;
        
        // Create alias in global scope
        let global_reg = q[0:2];
        
        // Create alias in local scope
        if ( true ) {
            let local_reg = q[3:4];
            h local_reg[0];  // Should apply H to q[3]
        }
        
        // Use global alias after local scope
        x global_reg[1];  // Should apply X to q[1]
        
        b = measure q[3];  // Measure q[3] which had H applied
        """

        # Create a simulator
        simulator = StateVectorSimulator(5, 1000)

        # Evolve the program using the branched simulator operators
        branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

        # Verify that we have two paths (from measuring q[3])
        @test length(branched_sim.active_paths) == 2
        
        # Calculate the current state for one path
        path_idx = branched_sim.active_paths[1]
        state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
        
        # Verify that X was applied to q[1]
        # Extract the probability for q[1] = 1
        prob_q1_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 3) != 0)
        
        @test isapprox(prob_q1_1, 1.0, atol=1e-10)
    end

    @testset "Error handling for physical qubits" begin
        # Create an OpenQASM program that tries to alias physical qubits
        qasm_source = """
        OPENQASM 3.0;
        
        // Try to alias a physical qubit
        let physical_alias = \$0;
        
        // This should cause an error
        h physical_alias;
        """

        # Create a simulator
        simulator = StateVectorSimulator(1, 1000)

        # Evolve the program should throw an error
        @test_throws ErrorException BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
    end

    @testset "Complex aliasing with measurements" begin
        # Create an OpenQASM program with aliases and measurements
        qasm_source = """
        OPENQASM 3.0;
        qubit[5] q;
        bit[2] b;
        
        // Create aliases
        let reg1 = q[0:1];
        let reg2 = q[2:3];
        
        // Apply gates using aliases
        h reg1[0];
        cnot reg1[0], reg2[0];
        
        // Measure using aliases
        b[0] = measure reg1[0];
        b[1] = measure reg2[0];
        
        // Conditional operations based on measurements
        if (b[0] == b[1]) {
            // Create a new alias after measurement
            let result_reg = q[4:4];
            x result_reg[0];  // Should apply X to q[4] if measurements match
        }
        """

        # Create a simulator
        simulator = StateVectorSimulator(5, 1000)

        # Evolve the program using the branched simulator operators
        branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

        # We should have 2 paths (since we created an entangled state, measurements should be correlated)
        @test length(branched_sim.active_paths) == 2
        
        # Find paths where measurements match
        matching_path = nothing
        
        for path_idx in branched_sim.active_paths
            measurements = branched_sim.measurements[path_idx]
            if measurements["reg1[0]"][1] == measurements["reg2[0]"][1]
                matching_path = path_idx
                break
            end
        end
        
        # Verify that we found a path where measurements match
        @test !isnothing(matching_path)
        
        # Calculate the state for the matching path
        state = BraketSimulator.calculate_current_state(branched_sim, matching_path)
        
        # Verify that X was applied to q[4] in the matching path
        # Extract the probability for q[4] = 1
        prob_q4_1 = sum(abs2.(state[i]) for i in 1:32 if (i-1) & (1 << 4) != 0)
        
        @test isapprox(prob_q4_1, 1.0, atol=1e-10)
    end
end
