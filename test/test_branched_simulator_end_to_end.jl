using Test, LinearAlgebra, Random
using BraketSimulator

@testset "Branched Simulator End-to-End Tests" begin
    @testset "Simple Mid-Circuit Measurement and Conditional" begin
        # Create a simple OpenQASM program with mid-circuit measurement
        qasm_source = """
        OPENQASM 3.0;
        bit b;
        qubit[2] q;
        
        h q[0];       // Put qubit 0 in superposition
        b = measure q[0];  // Measure qubit 0
        if (b == 1) {  // Conditional on measurement
            x q[1];    // Apply X to qubit 1
        }
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(2, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # Verify that we have two paths (one for each measurement outcome)
        @test length(branched_sim.active_paths) == 2
        
        # Get the path indices
        path1_idx = branched_sim.active_paths[1]
        path2_idx = branched_sim.active_paths[2]
        
        # Verify that one path measured 0 and the other measured 1
        @test (branched_sim.measurements[path1_idx][0] == 0 && branched_sim.measurements[path2_idx][0] == 1) ||
              (branched_sim.measurements[path1_idx][0] == 1 && branched_sim.measurements[path2_idx][0] == 0)
        
        # Find which path measured 1
        path_with_1 = branched_sim.measurements[path1_idx][0] == 1 ? path1_idx : path2_idx
        path_with_0 = branched_sim.measurements[path1_idx][0] == 0 ? path1_idx : path2_idx
        
        # Verify that the conditional X gate was applied in the path where b=1
        # For the path where b=1, qubit 1 should be in state |1⟩
        # For the path where b=0, qubit 1 should be in state |0⟩
        
        # Get the statevector for the path where b=1
        state_with_1 = branched_sim.states[path_with_1]
        
        # Get the statevector for the path where b=0
        state_with_0 = branched_sim.states[path_with_0]
        
        # Check if qubit 1 is in state |1⟩ for the path where b=1
        # Since the statevector is reduced after measurement, we need to check the remaining qubit
        @test abs(state_with_1[2]) ≈ 1.0 atol=1e-10
        
        # Check if qubit 1 is in state |0⟩ for the path where b=0
        @test abs(state_with_0[1]) ≈ 1.0 atol=1e-10
        
        # Verify that the weights are approximately equal (since H gives 50/50)
        @test branched_sim.weights[path_with_0] ≈ 0.5 atol=0.1
        @test branched_sim.weights[path_with_1] ≈ 0.5 atol=0.1
        
        # Allocate shots and verify distribution
        shot_allocation = BraketSimulator.allocate_shots(branched_sim)
        @test shot_allocation[path_with_0] ≈ 500 atol=50
        @test shot_allocation[path_with_1] ≈ 500 atol=50
    end
    
    @testset "Multiple Measurements and Branching Paths" begin
        # Create an OpenQASM program with multiple measurements
        qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[3] q;
        
        h q[0];       // Put qubit 0 in superposition
        h q[1];       // Put qubit 1 in superposition
        b[0] = measure q[0];  // Measure qubit 0
        b[1] = measure q[1];  // Measure qubit 1
        
        if (b[0] == 1 && b[1] == 1) {  // Both measured 1
            x q[2];    // Apply X to qubit 2
        } else if (b[0] == 1) {  // Only first qubit measured 1
            h q[2];    // Apply H to qubit 2
        } else if (b[1] == 1) {  // Only second qubit measured 1
            z q[2];    // Apply Z to qubit 2
        }
        // If both measured 0, do nothing to qubit 2
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(3, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # Verify that we have four paths (one for each combination of measurement outcomes)
        @test length(branched_sim.active_paths) == 4
        
        # Verify that the weights are approximately equal (since H gives 50/50 for each qubit)
        for path_idx in branched_sim.active_paths
            @test branched_sim.weights[path_idx] ≈ 0.25 atol=0.1
        end
        
        # Find paths for each measurement combination
        path_00 = nothing
        path_01 = nothing
        path_10 = nothing
        path_11 = nothing
        
        for path_idx in branched_sim.active_paths
            measurements = branched_sim.measurements[path_idx]
            if measurements[0] == 0 && measurements[1] == 0
                path_00 = path_idx
            elseif measurements[0] == 0 && measurements[1] == 1
                path_01 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 0
                path_10 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 1
                path_11 = path_idx
            end
        end
        
        # Verify that we found all four paths
        @test !isnothing(path_00)
        @test !isnothing(path_01)
        @test !isnothing(path_10)
        @test !isnothing(path_11)
        
        # Verify the state of qubit 2 for each path
        
        # Path 00: No operation, should be |0⟩
        @test abs(branched_sim.states[path_00][1]) ≈ 1.0 atol=1e-10
        
        # Path 01: Z operation, should be |0⟩ (Z doesn't change |0⟩)
        @test abs(branched_sim.states[path_01][1]) ≈ 1.0 atol=1e-10
        
        # Path 10: H operation, should be (|0⟩ + |1⟩)/√2
        @test abs(branched_sim.states[path_10][1]) ≈ 1/sqrt(2) atol=1e-10
        @test abs(branched_sim.states[path_10][2]) ≈ 1/sqrt(2) atol=1e-10
        
        # Path 11: X operation, should be |1⟩
        @test abs(branched_sim.states[path_11][2]) ≈ 1.0 atol=1e-10
        
        # Allocate shots and verify distribution
        shot_allocation = BraketSimulator.allocate_shots(branched_sim)
        @test shot_allocation[path_00] ≈ 250 atol=50
        @test shot_allocation[path_01] ≈ 250 atol=50
        @test shot_allocation[path_10] ≈ 250 atol=50
        @test shot_allocation[path_11] ≈ 250 atol=50
    end
    
    @testset "Complex Conditional Logic" begin
        # Create an OpenQASM program with complex conditional logic
        qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[3] q;
        
        h q[0];       // Put qubit 0 in superposition
        h q[1];       // Put qubit 1 in superposition
        
        b[0] = measure q[0];  // Measure qubit 0
        
        if (b[0] == 0) {
            h q[1];    // Apply H to qubit 1 if qubit 0 measured 0
        }
        
        b[1] = measure q[1];  // Measure qubit 1
        
        // Nested conditionals
        if (b[0] == 1) {
            if (b[1] == 1) {
                x q[2];    // Apply X to qubit 2 if both measured 1
            } else {
                h q[2];    // Apply H to qubit 2 if q0=1, q1=0
            }
        } else {
            if (b[1] == 1) {
                z q[2];    // Apply Z to qubit 2 if q0=0, q1=1
            } else {
                // Do nothing if both measured 0
            }
        }
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(3, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # Verify that we have four paths (one for each combination of measurement outcomes)
        @test length(branched_sim.active_paths) == 4
        
        # Find paths for each measurement combination
        path_00 = nothing
        path_01 = nothing
        path_10 = nothing
        path_11 = nothing
        
        for path_idx in branched_sim.active_paths
            measurements = branched_sim.measurements[path_idx]
            if measurements[0] == 0 && measurements[1] == 0
                path_00 = path_idx
            elseif measurements[0] == 0 && measurements[1] == 1
                path_01 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 0
                path_10 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 1
                path_11 = path_idx
            end
        end
        
        # Verify that we found all four paths
        @test !isnothing(path_00)
        @test !isnothing(path_01)
        @test !isnothing(path_10)
        @test !isnothing(path_11)
        
        # Verify the state of qubit 2 for each path
        
        # Path 00: No operation, should be |0⟩
        @test abs(branched_sim.states[path_00][1]) ≈ 1.0 atol=1e-10
        
        # Path 01: Z operation, should be |0⟩ (Z doesn't change |0⟩)
        @test abs(branched_sim.states[path_01][1]) ≈ 1.0 atol=1e-10
        
        # Path 10: H operation, should be (|0⟩ + |1⟩)/√2
        @test abs(branched_sim.states[path_10][1]) ≈ 1/sqrt(2) atol=1e-10
        @test abs(branched_sim.states[path_10][2]) ≈ 1/sqrt(2) atol=1e-10
        
        # Path 11: X operation, should be |1⟩
        @test abs(branched_sim.states[path_11][2]) ≈ 1.0 atol=1e-10
    end
    
    @testset "Classical Variable Manipulation" begin
        # Create an OpenQASM program with classical variable manipulation
        qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[3] q;
        int[32] count = 0;
        
        h q[0];       // Put qubit 0 in superposition
        h q[1];       // Put qubit 1 in superposition
        
        b[0] = measure q[0];  // Measure qubit 0
        b[1] = measure q[1];  // Measure qubit 1
        
        // Update count based on measurements
        if (b[0] == 1) {
            count = count + 1;
        }
        if (b[1] == 1) {
            count = count + 1;
        }
        
        // Apply operations based on count
        if (count == 0) {
            // Do nothing
        } else if (count == 1) {
            h q[2];    // Apply H to qubit 2 if one qubit measured 1
        } else if (count == 2) {
            x q[2];    // Apply X to qubit 2 if both qubits measured 1
        }
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(3, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # Verify that we have four paths (one for each combination of measurement outcomes)
        @test length(branched_sim.active_paths) == 4
        
        # Find paths for each measurement combination
        path_00 = nothing
        path_01 = nothing
        path_10 = nothing
        path_11 = nothing
        
        for path_idx in branched_sim.active_paths
            measurements = branched_sim.measurements[path_idx]
            if measurements[0] == 0 && measurements[1] == 0
                path_00 = path_idx
            elseif measurements[0] == 0 && measurements[1] == 1
                path_01 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 0
                path_10 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 1
                path_11 = path_idx
            end
        end
        
        # Verify that we found all four paths
        @test !isnothing(path_00)
        @test !isnothing(path_01)
        @test !isnothing(path_10)
        @test !isnothing(path_11)
        
        # Verify the count variable for each path
        @test BraketSimulator.get_variable(branched_sim, path_00, "count") == 0
        @test BraketSimulator.get_variable(branched_sim, path_01, "count") == 1
        @test BraketSimulator.get_variable(branched_sim, path_10, "count") == 1
        @test BraketSimulator.get_variable(branched_sim, path_11, "count") == 2
        
        # Verify the state of qubit 2 for each path
        
        # Path 00: count=0, no operation, should be |0⟩
        @test abs(branched_sim.states[path_00][1]) ≈ 1.0 atol=1e-10
        
        # Path 01: count=1, H operation, should be (|0⟩ + |1⟩)/√2
        @test abs(branched_sim.states[path_01][1]) ≈ 1/sqrt(2) atol=1e-10
        @test abs(branched_sim.states[path_01][2]) ≈ 1/sqrt(2) atol=1e-10
        
        # Path 10: count=1, H operation, should be (|0⟩ + |1⟩)/√2
        @test abs(branched_sim.states[path_10][1]) ≈ 1/sqrt(2) atol=1e-10
        @test abs(branched_sim.states[path_10][2]) ≈ 1/sqrt(2) atol=1e-10
        
        # Path 11: count=2, X operation, should be |1⟩
        @test abs(branched_sim.states[path_11][2]) ≈ 1.0 atol=1e-10
    end
    
    @testset "Loop Dependent on Measurement Results" begin
        # Create an OpenQASM program with a loop dependent on measurement results
        qasm_source = """
        OPENQASM 3.0;
        bit b;
        qubit[2] q;
        int[32] count = 0;
        
        // Initialize qubit 0 to |0⟩
        reset q[0];
        
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
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(2, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # We should have 4 paths:
        # 1. Measured 1 on first attempt
        # 2. Measured 0 on first attempt, then 1 on second attempt
        # 3. Measured 0 on first and second attempts, then 1 on third attempt
        # 4. Measured 0 on all three attempts
        @test length(branched_sim.active_paths) == 4
        
        # Find paths for each count value
        paths_by_count = Dict{Int, Vector{Int}}()
        
        for path_idx in branched_sim.active_paths
            count = BraketSimulator.get_variable(branched_sim, path_idx, "count")
            if !haskey(paths_by_count, count)
                paths_by_count[count] = Int[]
            end
            push!(paths_by_count[count], path_idx)
        end
        
        # Verify that we have paths for counts 1, 2, 3
        @test haskey(paths_by_count, 1)
        @test haskey(paths_by_count, 2)
        @test haskey(paths_by_count, 3)
        
        # For paths with count < 3, b should be 1 and qubit 1 should be in state |1⟩
        for count in [1, 2]
            for path_idx in paths_by_count[count]
                b = BraketSimulator.get_variable(branched_sim, path_idx, "b")
                @test b == 1
                @test abs(branched_sim.states[path_idx][2]) ≈ 1.0 atol=1e-10
            end
        end
        
        # For paths with count = 3, check each path individually
        for path_idx in paths_by_count[3]
            b = BraketSimulator.get_variable(branched_sim, path_idx, "b")
            if b == 1
                # If b=1, qubit 1 should be in state |1⟩
                @test abs(branched_sim.states[path_idx][2]) ≈ 1.0 atol=1e-10
            else
                # If b=0, qubit 1 should be in state |0⟩
                @test abs(branched_sim.states[path_idx][1]) ≈ 1.0 atol=1e-10
            end
        end
        
        # Verify the weights add up to approximately 1
        total_weight = sum(branched_sim.weights[path_idx] for path_idx in branched_sim.active_paths)
        @test total_weight ≈ 1.0 atol=1e-10
    end
    
    @testset "Memory Optimization with Large Qubit Count" begin
        # Create an OpenQASM program with a large number of qubits
        # We'll measure half of them to test memory optimization
        n_qubits = 10  # 10 qubits = 2^10 = 1024 amplitudes
        half_qubits = div(n_qubits, 2)
        
        # Build the QASM string programmatically
        qasm_source = """
        OPENQASM 3.0;
        bit[$half_qubits] b;
        qubit[$n_qubits] q;
        
        // Put all qubits in superposition
        h q;
        
        // Measure half of the qubits
        """
        
        # Add measurement operations for half of the qubits
        for i in 0:(half_qubits-1)
            qasm_source *= "b[$i] = measure q[$i];\n"
        end
        
        # Add a conditional operation based on the first measurement
        qasm_source *= """
        
        // Apply X to the last qubit if the first measurement is 1
        if (b[0] == 1) {
            x q[$n_qubits-1];
        }
        """
        
        # Create a simulator with 100 shots (fewer shots for large qubit count)
        simulator = StateVectorSimulator(n_qubits, 100)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # We should have 2^half_qubits paths (one for each combination of measurement outcomes)
        @test length(branched_sim.active_paths) == 2^half_qubits
        
        # Verify that the statevector size is reduced
        for path_idx in branched_sim.active_paths
            # After measuring half_qubits qubits, the statevector should have 2^(n_qubits - half_qubits) amplitudes
            @test length(branched_sim.states[path_idx]) == 2^(n_qubits - half_qubits)
            
            # Verify that the measured qubits are recorded as removed
            for i in 0:(half_qubits-1)
                @test haskey(branched_sim.removed_qubits[path_idx], i)
            end
            
            # Verify that the qubit mapping is updated
            for i in 0:(half_qubits-1)
                @test !haskey(branched_sim.qubit_mapping[path_idx], i)
            end
            for i in half_qubits:(n_qubits-1)
                @test haskey(branched_sim.qubit_mapping[path_idx], i)
            end
        end
        
        # Find paths where the first measurement is 0 or 1
        paths_with_0 = Int[]
        paths_with_1 = Int[]
        
        for path_idx in branched_sim.active_paths
            if branched_sim.measurements[path_idx][0] == 0
                push!(paths_with_0, path_idx)
            else
                push!(paths_with_1, path_idx)
            end
        end
        
        # Verify that we have paths for both measurement outcomes
        @test !isempty(paths_with_0)
        @test !isempty(paths_with_1)
        
        # For paths where the first measurement is 1, the last qubit should be in state |1⟩
        for path_idx in paths_with_1
            # Get the mapped index of the last qubit
            last_qubit_mapped = branched_sim.qubit_mapping[path_idx][n_qubits-1]
            
            # Calculate the index in the reduced statevector where the last qubit is 1
            # This is more complex due to the reduced statevector
            one_state_index = 1 << (n_qubits - half_qubits - 1 - last_qubit_mapped)
            
            # The amplitude at this index should be 1
            # We need to check all basis states where the last qubit is 1
            total_prob = 0.0
            for i in 0:(2^(n_qubits - half_qubits)-1)
                if (i & one_state_index) != 0  # Check if the bit corresponding to the last qubit is 1
                    total_prob += abs2(branched_sim.states[path_idx][i+1])
                end
            end
            @test total_prob ≈ 1.0 atol=1e-10
        end
        
        # For paths where the first measurement is 0, the last qubit should be in state |0⟩
        for path_idx in paths_with_0
            # Get the mapped index of the last qubit
            last_qubit_mapped = branched_sim.qubit_mapping[path_idx][n_qubits-1]
            
            # Calculate the index in the reduced statevector where the last qubit is 0
            # This is more complex due to the reduced statevector
            zero_state_index = 1 << (n_qubits - half_qubits - 1 - last_qubit_mapped)
            
            # The amplitude at this index should be 1
            # We need to check all basis states where the last qubit is 0
            total_prob = 0.0
            for i in 0:(2^(n_qubits - half_qubits)-1)
                if (i & zero_state_index) == 0  # Check if the bit corresponding to the last qubit is 0
                    total_prob += abs2(branched_sim.states[path_idx][i+1])
                end
            end
            @test total_prob ≈ 1.0 atol=1e-10
        end
    end
    
    @testset "Multiple Measurements on Same Qubit" begin
        # Create an OpenQASM program with multiple measurements on the same qubit
        qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[2] q;
        
        // Put qubit 0 in superposition
        h q[0];
        
        // First measurement
        b[0] = measure q[0];
        
        // Apply X to qubit 0 if measured 0
        if (b[0] == 0) {
            x q[0];
        }
        
        // Second measurement (should always be 1)
        b[1] = measure q[0];
        
        // Apply X to qubit 1 if both measurements are the same
        if (b[0] == b[1]) {
            x q[1];
        }
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(2, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # We should have 2 paths (one for each outcome of the first measurement)
        @test length(branched_sim.active_paths) == 2
        
        # Find paths for each first measurement outcome
        path_0 = nothing
        path_1 = nothing
        
        for path_idx in branched_sim.active_paths
            if branched_sim.measurements[path_idx][0] == 0
                path_0 = path_idx
            else
                path_1 = path_idx
            end
        end
        
        # Verify that we found both paths
        @test !isnothing(path_0)
        @test !isnothing(path_1)
        
        # Verify that the second measurement is always 1
        @test branched_sim.measurements[path_0][0] == 1
        @test branched_sim.measurements[path_1][0] == 1
        
        # For path_0, b[0]=0 and b[1]=1, so they're different, and qubit 1 should be |0⟩
        @test abs(branched_sim.states[path_0][1]) ≈ 1.0 atol=1e-10
        
        # For path_1, b[0]=1 and b[1]=1, so they're the same, and qubit 1 should be |1⟩
        @test abs(branched_sim.states[path_1][2]) ≈ 1.0 atol=1e-10
    end
    
    @testset "Quantum Teleportation" begin
        # Create an OpenQASM program for quantum teleportation
        qasm_source = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[3] q;
        
        // Prepare the state to teleport on qubit 0
        // Let's use |+⟩ state
        h q[0];
        
        // Create Bell pair between qubits 1 and 2
        h q[1];
        cx q[1], q[2];
        
        // Perform teleportation protocol
        cx q[0], q[1];
        h q[0];
        b[0] = measure q[0];
        b[1] = measure q[1];
        
        // Apply corrections based on measurement results
        if (b[1] == 1) {
            x q[2];  // Apply Pauli X
        }
        if (b[0] == 1) {
            z q[2];  // Apply Pauli Z
        }
        
        // At this point, qubit 2 should be in the |+⟩ state
        """
        
        # Create a simulator with 1000 shots
        simulator = StateVectorSimulator(3, 1000)
        
        # Evolve the program using the branched simulator
        branched_sim = BraketSimulator.evolve_branched!(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
        
        # We should have 4 paths (one for each combination of measurement outcomes)
        @test length(branched_sim.active_paths) == 4
        
        # Find paths for each measurement combination
        path_00 = nothing
        path_01 = nothing
        path_10 = nothing
        path_11 = nothing
        
        for path_idx in branched_sim.active_paths
            measurements = branched_sim.measurements[path_idx]
            if measurements[0] == 0 && measurements[1] == 0
                path_00 = path_idx
            elseif measurements[0] == 0 && measurements[1] == 1
                path_01 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 0
                path_10 = path_idx
            elseif measurements[0] == 1 && measurements[1] == 1
                path_11 = path_idx
            end
        end
        
        # Verify that we found all four paths
        @test !isnothing(path_00)
        @test !isnothing(path_01)
        @test !isnothing(path_10)
        @test !isnothing(path_11)
        
        # For all paths, qubit 2 should be in the |+⟩ state (|0⟩ + |1⟩)/√2
        # regardless of the measurement outcomes, due to the corrections
        for path_idx in [path_00, path_01, path_10, path_11]
            # Get the statevector for this path
            state = branched_sim.states[path_idx]
            
            # Check if qubit 2 is in state |+⟩
            @test abs(state[1]) ≈ 1/sqrt(2) atol=1e-10
            @test abs(state[2]) ≈ 1/sqrt(2) atol=1e-10
        end
    end
end
