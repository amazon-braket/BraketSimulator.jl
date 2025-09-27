using Test, LinearAlgebra, Random
using BraketSimulator

@testset "Reset Operator Tests" begin
    @testset "1. Basic Reset Operation" begin
        @testset "1.1 Reset a single qubit from |0⟩ state" begin
            # Create a simple OpenQASM program with reset on |0⟩ state
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            reset q[0];  // Reset qubit 0 (already in |0⟩ state)
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that the qubit is in |0⟩ state
            @test abs(state[1]) ≈ 1.0 atol=1e-10
            @test abs(state[2]) ≈ 0.0 atol=1e-10
        end

        @testset "1.2 Reset a single qubit from |1⟩ state" begin
            # Create a simple OpenQASM program with reset on |1⟩ state
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            x q[0];      // Put qubit 0 in |1⟩ state
            reset q[0];  // Reset qubit 0 to |0⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that the qubit is in |0⟩ state after reset
            @test abs(state[1]) ≈ 1.0 atol=1e-10
            @test abs(state[2]) ≈ 0.0 atol=1e-10
        end

        @testset "1.3 Reset a single qubit from superposition" begin
            # Create a simple OpenQASM program with reset on superposition state
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            h q[0];      // Put qubit 0 in superposition
            reset q[0];  // Reset qubit 0 to |0⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that the qubit is in |0⟩ state after reset
            @test abs(state[1]) ≈ 1.0 atol=1e-10
            @test abs(state[2]) ≈ 0.0 atol=1e-10
        end
    end

    @testset "2. Reset in Multi-Qubit Systems" begin
        @testset "2.1 Reset one qubit in a two-qubit system" begin
            # Create an OpenQASM program with reset on one qubit in a two-qubit system
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            x q[0];      // Put qubit 0 in |1⟩ state
            x q[1];      // Put qubit 1 in |1⟩ state
            reset q[0];  // Reset qubit 0 to |0⟩ state, qubit 1 should remain in |1⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that qubit 0 is in |0⟩ state and qubit 1 is in |1⟩ state
            # The state should be |01⟩
            @test abs(state[1]) ≈ 0.0 atol=1e-10
            @test abs(state[2]) ≈ 1.0 atol=1e-10
            @test abs(state[3]) ≈ 0.0 atol=1e-10
            @test abs(state[4]) ≈ 0.0 atol=1e-10
        end

        @testset "2.2 Reset both qubits in a two-qubit system" begin
            # Create an OpenQASM program with reset on both qubits in a two-qubit system
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            x q[0];      // Put qubit 0 in |1⟩ state
            x q[1];      // Put qubit 1 in |1⟩ state
            reset q[0];  // Reset qubit 0 to |0⟩ state
            reset q[1];  // Reset qubit 1 to |0⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that both qubits are in |0⟩ state
            # The state should be |00⟩
            @test abs(state[1]) ≈ 1.0 atol=1e-10
            @test abs(state[2]) ≈ 0.0 atol=1e-10
            @test abs(state[3]) ≈ 0.0 atol=1e-10
            @test abs(state[4]) ≈ 0.0 atol=1e-10
        end

        @testset "2.3 Reset in an entangled system" begin
            # Create an OpenQASM program with reset on an entangled system
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            h q[0];          // Put qubit 0 in superposition
            cnot q[0], q[1]; // Entangle qubits 0 and 1
            reset q[0];      // Reset qubit 0 to |0⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # After reset, the entanglement should be broken
            # The state should be a mixture of |00⟩ and |01⟩ with equal probability
            # But since we're using a state vector simulator, we should get a pure state
            # The exact state depends on the implementation details of the reset operation
            # We'll check that qubit 0 is definitely in |0⟩ state
            @test abs(state[1])^2 + abs(state[2])^2 ≈ 1.0 atol=1e-10
            @test abs(state[3]) ≈ 0.0 atol=1e-10
            @test abs(state[4]) ≈ 0.0 atol=1e-10
        end
    end

    @testset "3. Reset with Subsequent Operations" begin
        @testset "3.1 Operations after reset" begin
            # Create an OpenQASM program with operations after reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            x q[0];      // Put qubit 0 in |1⟩ state
            x q[1];      // Put qubit 1 in |1⟩ state
            reset q[0];  // Reset qubit 0 to |0⟩ state
            h q[0];      // Apply H to qubit 0
            cnot q[0], q[1]; // Apply CNOT with qubit 0 as control
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # After reset, H, and CNOT, the state should be (|01⟩ + |10⟩)/√2
            @test abs(state[1]) ≈ 0.0 atol=1e-10
            @test abs(state[2]) ≈ 1/sqrt(2) atol=1e-10
            @test abs(state[3]) ≈ 1/sqrt(2) atol=1e-10
            @test abs(state[4]) ≈ 0.0 atol=1e-10
        end

        @testset "3.2 Reset between operations" begin
            # Create an OpenQASM program with reset between operations
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            h q[0];      // Put qubit 0 in superposition
            cnot q[0], q[1]; // Entangle qubits 0 and 1
            reset q[0];  // Reset qubit 0 to |0⟩ state
            h q[0];      // Apply H to qubit 0
            cnot q[0], q[1]; // Apply CNOT with qubit 0 as control
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # After the second CNOT, the state should be a superposition
            # The exact state depends on the implementation details of the reset operation
            # We'll check that the state is normalized
            @test abs(state[1])^2 ≈ 0.25 atol=1e-10
            @test abs(state[2])^2 ≈ 0.25 atol=1e-10
            @test abs(state[3])^2 ≈ 0.25 atol=1e-10
            @test abs(state[4])^2 ≈ 0.25 atol=1e-10
        end
    end

    @testset "4. Reset with Measurements" begin
        @testset "4.1 Measurement after reset" begin
            # Create an OpenQASM program with measurement after reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            bit[1] b;
            x q[0];      // Put qubit 0 in |1⟩ state
            reset q[0];  // Reset qubit 0 to |0⟩ state
            b[0] = measure q[0]; // Measure qubit 0
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that we have only one path (since the measurement outcome is deterministic)
            @test length(branched_sim.active_paths) == 1

            # Verify that the measurement outcome is 0
            path_idx = branched_sim.active_paths[1]
            @test branched_sim.measurements[path_idx]["q[0]"][1] == 0
        end

        @testset "4.2 Reset after measurement" begin
            # Create an OpenQASM program with reset after measurement
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            bit[1] b;
            h q[0];      // Put qubit 0 in superposition
            b[0] = measure q[0]; // Measure qubit 0
            reset q[0];  // Reset qubit 0 to |0⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that we have two paths (one for each measurement outcome)
            @test length(branched_sim.active_paths) == 2

            # For each path, verify that the final state is |0⟩
            for path_idx in branched_sim.active_paths
                state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
                @test abs(state[1]) ≈ 1.0 atol=1e-10
                @test abs(state[2]) ≈ 0.0 atol=1e-10
            end
        end

        @testset "4.3 Reset and re-measure" begin
            # Create an OpenQASM program with reset and re-measure
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            bit[2] b;
            h q[0];          // Put qubit 0 in superposition
            b[0] = measure q[0]; // Measure qubit 0
            reset q[0];      // Reset qubit 0 to |0⟩ state
            b[1] = measure q[0]; // Measure qubit 0 again
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that we have two paths (one for each outcome of the first measurement)
            @test length(branched_sim.active_paths) == 2

            # For each path, verify that the second measurement outcome is 0
            for path_idx in branched_sim.active_paths
                @test branched_sim.measurements[path_idx]["q[0]"][2] == 0
            end
        end
    end

    @testset "5. Reset in Conditional Blocks" begin
        @testset "5.1 Conditional reset based on measurement" begin
            # Create an OpenQASM program with conditional reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            bit[2] b;
            h q[0];          // Put qubit 0 in superposition
            x q[1];          // Put qubit 1 in |1⟩ state
            b[0] = measure q[0]; // Measure qubit 0
            if (b[0] == 1) {
                reset q[1];  // Reset qubit 1 to |0⟩ state if qubit 0 was measured as 1
            }
            b[1] = measure q[1]; // Measure qubit 1
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that we have two paths (one for each outcome of the first measurement)
            @test length(branched_sim.active_paths) == 2

            # Find paths for each measurement outcome
            path_0 = nothing
            path_1 = nothing

            for path_idx in branched_sim.active_paths
                if branched_sim.measurements[path_idx]["q[0]"][1] == 0
                    path_0 = path_idx
                else
                    path_1 = path_idx
                end
            end

            # Verify that we found both paths
            @test !isnothing(path_0)
            @test !isnothing(path_1)

            # For path_0 (b[0] = 0), qubit 1 should remain in |1⟩ state
            @test branched_sim.measurements[path_0]["q[1]"][1] == 1

            # For path_1 (b[0] = 1), qubit 1 should be reset to |0⟩ state
            @test branched_sim.measurements[path_1]["q[1]"][1] == 0
        end
    end

    @testset "6. Reset in Complex Circuits" begin
        @testset "6.1 Reset in quantum teleportation" begin
            # Create an OpenQASM program for quantum teleportation with reset
            qasm_source = """
            OPENQASM 3.0;
            bit[2] b;
            qubit[3] q;

            // Prepare the state to teleport on qubit 0
            // Let's use |+⟩ state
            h q[0];

            // Create Bell pair between qubits 1 and 2
            h q[1];
            cnot q[1], q[2];

            // Perform teleportation protocol
            cnot q[0], q[1];
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

            // Reset qubits 0 and 1 (they are no longer needed)
            reset q[0];
            reset q[1];

            // At this point, qubit 2 should be in the |+⟩ state
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(3, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that we have 4 paths (one for each combination of measurement outcomes)
            @test length(branched_sim.active_paths) == 4

            # For each path, verify that qubits 0 and 1 are in |0⟩ state and qubit 2 is in |+⟩ state
            for path_idx in branched_sim.active_paths
                state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

                # Extract the reduced state of qubit 2
                # Since qubits 0 and 1 are in |0⟩ state, we only need to check state[1] and state[2]
                @test abs(state[1]) ≈ 1/sqrt(2) atol=1e-10
                @test abs(state[2]) ≈ 1/sqrt(2) atol=1e-10
                @test abs(state[3]) ≈ 0.0 atol=1e-10
                @test abs(state[4]) ≈ 0.0 atol=1e-10
                @test abs(state[5]) ≈ 0.0 atol=1e-10
                @test abs(state[6]) ≈ 0.0 atol=1e-10
                @test abs(state[7]) ≈ 0.0 atol=1e-10
                @test abs(state[8]) ≈ 0.0 atol=1e-10
            end
        end

        @testset "6.2 Reset in error correction" begin
            # Create an OpenQASM program for a simple bit-flip code with reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[4] q;  // 3 data qubits + 1 ancilla
            bit[3] syndrome;
            bit data;

            // Prepare a state to protect
            h q[0];

            // Encode using a simple repetition code
            cnot q[0], q[1];
            cnot q[0], q[2];

            // Introduce an error (bit flip on qubit 1)
            x q[1];

            // Detect the error using an ancilla
            reset q[3];  // Reset ancilla to |0⟩
            cnot q[0], q[3];
            cnot q[1], q[3];
            syndrome[0] = measure q[3];

            reset q[3];  // Reset ancilla for next syndrome
            cnot q[1], q[3];
            cnot q[2], q[3];
            syndrome[1] = measure q[3];

            // Correct the error based on syndrome
            if (syndrome[0] == 1 && syndrome[1] == 1) {
                x q[1];  // Correct bit flip on qubit 1
            }
            else if (syndrome[0] == 1 && syndrome[1] == 0) {
                x q[0];  // Correct bit flip on qubit 0
            }
            else if (syndrome[0] == 0 && syndrome[1] == 1) {
                x q[2];  // Correct bit flip on qubit 2
            }

            // Decode
            cnot q[0], q[1];
            cnot q[0], q[2];

            // Measure the result
            data = measure q[0];
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(4, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Verify that the error was corrected
            # The final measurement should be in a superposition state
            # We should have 2 paths (one for each possible measurement outcome of q[0])
            @test length(branched_sim.active_paths) == 2

            # Count the number of paths with each measurement outcome
            count_0 = 0
            count_1 = 0

            for path_idx in branched_sim.active_paths
                if branched_sim.measurements[path_idx]["q[0]"][1] == 0
                    count_0 += 1
                else
                    count_1 += 1
                end
            end

            # Both outcomes should be possible
            @test count_0 > 0
            @test count_1 > 0
        end
    end

    @testset "7. Reset Implementation Details" begin
        @testset "7.1 Reset preserves normalization" begin
            # Create an OpenQASM program to test normalization after reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[2] q;
            
            // Create a superposition state
            h q[0];
            h q[1];
            
            // Reset one qubit
            reset q[0];
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(2, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that the state is normalized
            @test sum(abs2.(state)) ≈ 1.0 atol=1e-10
        end

        @testset "7.2 Reset transfers amplitudes correctly" begin
            # Create an OpenQASM program to test amplitude transfer during reset
            qasm_source = """
            OPENQASM 3.0;
            qubit[1] q;
            
            // Create a specific state with known amplitudes
            ry(π/3) q[0];  // Creates cos(π/6)|0⟩ + sin(π/6)|1⟩
            
            // Reset the qubit
            reset q[0];
            """

            # Create a simulator
            simulator = BranchedSimulatorOperators(StateVectorSimulator(1, 1000))

            # Evolve the program using the branched simulator operators
            branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

            # Calculate the final state
            state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

            # Verify that the qubit is in |0⟩ state
            @test abs(state[1]) ≈ 1.0 atol=1e-10
            @test abs(state[2]) ≈ 0.0 atol=1e-10
        end
    end
end
