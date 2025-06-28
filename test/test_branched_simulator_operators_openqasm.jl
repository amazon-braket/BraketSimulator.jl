using Test, LinearAlgebra, Random
using BraketSimulator

@testset "Branched Simulator Operators with OpenQASM" begin
	@testset "1. Basic OpenQASM 3 Features" begin
		@testset "1.1 Basic initialization and simple operations" begin
			# Create a simple OpenQASM program with basic operations
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;

			h q[0];       // Put qubit 0 in superposition
			cnot q[0], q[1];  // Create Bell state
			"""

			# Parse the program
			circuit = BraketSimulator.to_circuit(qasm_source)

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Create a branched simulator operators instance
			branched_sim = BranchedSimulatorOperators(simulator)

			# Verify initial state
			@test length(branched_sim.instruction_sequences) == 1
			@test length(branched_sim.active_paths) == 1
			@test branched_sim.n_qubits == 0
			@test isempty(branched_sim.measurements[1])

			# Evolve the circuit
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that the instruction sequence contains the expected operations
			@test length(branched_sim.instruction_sequences[1]) == 2
		end

		@testset "1.2 Empty Circuit" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			@test length(branched_sim.active_paths) == 1
			@test isempty(branched_sim.instruction_sequences[1])
		end
	end

	@testset "2. Measurement Operations" begin
		@testset "2.1 Mid-circuit measurement" begin
			# Create an OpenQASM program with mid-circuit measurement
			qasm_source = """
			OPENQASM 3.0;
			bit b;
			qubit[2] q;

			h q[0];       // Put qubit 0 in superposition
			b = measure q[0];  // Measure qubit 0
			"""

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2

			# Get the path indices
			path1_idx = branched_sim.active_paths[1]
			path2_idx = branched_sim.active_paths[2]

			# Verify that one path measured 0 and the other measured 1
			@test (branched_sim.measurements[path1_idx]["q[0]"][1] == 0 && branched_sim.measurements[path2_idx]["q[0]"][1] == 1) ||
				  (branched_sim.measurements[path1_idx]["q[0]"][1] == 1 && branched_sim.measurements[path2_idx]["q[0]"][1] == 0)
		end

		@testset "2.2 Multiple measurements on same qubit" begin
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

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 2 paths (one for each outcome of the first measurement)
			@test length(branched_sim.active_paths) == 2

			# Find paths for each first measurement outcome
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

			# Verify that the second measurement is always 1
			@test branched_sim.measurements[path_0]["q[0]"][2] == 1
			@test branched_sim.measurements[path_1]["q[0]"][2] == 1

			# Calculate the final states for both paths
			state_0 = BraketSimulator.calculate_current_state(branched_sim, path_0)
			state_1 = BraketSimulator.calculate_current_state(branched_sim, path_1)

			# For path_0, b[0]=0 and b[1]=1, so they're different, and qubit 1 should be |0⟩
			@test abs(state_0[3]) ≈ 1.0 atol=1e-10

			# For path_1, b[0]=1 and b[1]=1, so they're the same, and qubit 1 should be |1⟩
			@test abs(state_1[4]) ≈ 1.0 atol=1e-10
		end
	end

	@testset "3. Conditional Operations" begin
		@testset "3.1 Simple conditional operations (feedforward)" begin
			# Create an OpenQASM program with conditional operations
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

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2

			# Find which path measured 1
			path_with_1 = nothing
			path_with_0 = nothing

			for path_idx in branched_sim.active_paths
				if branched_sim.measurements[path_idx]["q[0]"][1] == 1
					path_with_1 = path_idx
				else
					path_with_0 = path_idx
				end
			end

			@test !isnothing(path_with_1)
			@test !isnothing(path_with_0)

			# Verify that the X gate was applied in the path where b=1
			# The path with b=1 should have one more instruction (the X gate)
			@test length(branched_sim.instruction_sequences[path_with_1]) > length(branched_sim.instruction_sequences[path_with_0])

			# Calculate the final states for both paths
			state_with_1 = BraketSimulator.calculate_current_state(branched_sim, path_with_1)
			state_with_0 = BraketSimulator.calculate_current_state(branched_sim, path_with_0)

			# For the path where b=1, qubit 1 should be in state |1⟩
			# For the path where b=0, qubit 1 should be in state |0⟩
			@test abs(state_with_1[4]) ≈ 1.0 atol=1e-10  # |11⟩ state
			@test abs(state_with_0[1]) ≈ 1.0 atol=1e-10  # |00⟩ state
		end

		@testset "3.2 Complex conditional logic" begin
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

			# Create a simulator
			simulator = StateVectorSimulator(3, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have three paths
			@test length(branched_sim.active_paths) == 3

			# Find paths for each measurement combination
			path_00 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all paths
			@test !isnothing(path_00)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# Verify the state of qubit 2 for each path
			# Path 00: No operation, should be |0⟩
			@test abs(state_00[1]) ≈ 1.0 atol=1e-10

			# Path 10: H operation, should be (|0⟩ + |1⟩)/√2
			@test abs(state_10[5]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_10[6]) ≈ 1/sqrt(2) atol=1e-10

			# Path 11: X operation, should be |1⟩
			@test abs(state_11[8]) ≈ 1.0 atol=1e-10
		end

		@testset "3.3 Multiple measurements and branching paths" begin
			# Create an OpenQASM program with multiple measurements
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
				} else {
					h q[2];    // Apply H to qubit 2
				}
			} else {
				if (b[1] == 1) {  // Only second qubit measured 1
					z q[2];    // Apply Z to qubit 2
				}
			}
			// If both measured 0, do nothing to qubit 2
			"""

			# Create a simulator
			simulator = StateVectorSimulator(3, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have four paths (one for each combination of measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# Find paths for each measurement combination
			path_00 = nothing
			path_01 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 1
					path_01 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all four paths
			@test !isnothing(path_00)
			@test !isnothing(path_01)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_01 = BraketSimulator.calculate_current_state(branched_sim, path_01)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# Verify the state of qubit 2 for each path
			# Path 00: No operation, should be |0⟩
			@test abs(state_00[1]) ≈ 1.0 atol=1e-10

			# Path 01: Z operation, should be |0⟩ (Z doesn't change |0⟩)
			@test abs(state_01[3]) ≈ 1.0 atol=1e-10

			# Path 10: H operation, should be (|0⟩ + |1⟩)/√2
			@test abs(state_10[5]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_10[6]) ≈ 1/sqrt(2) atol=1e-10

			# Path 11: X operation, should be |1⟩
			@test abs(state_11[8]) ≈ 1.0 atol=1e-10
		end
	end

	@testset "4. Classical Data Types and Operations" begin
		@testset "4.1 Classical variable manipulation" begin
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
			switch (count) {
			case 1 {
				h q[2];    // Apply H to qubit 2 if one qubit measured 1
			}
			case 2 {
				x q[2];    // Apply X to qubit 2 if both qubits measured 1
			}
			}
			"""

			# Create a simulator
			simulator = StateVectorSimulator(3, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have four paths (one for each combination of measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# Find paths for each measurement combination
			path_00 = nothing
			path_01 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 1
					path_01 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all four paths
			@test !isnothing(path_00)
			@test !isnothing(path_01)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Verify the count variable for each path
			@test BraketSimulator.get_variable(branched_sim, path_00, "count").val == 0
			@test BraketSimulator.get_variable(branched_sim, path_01, "count").val == 1
			@test BraketSimulator.get_variable(branched_sim, path_10, "count").val == 1
			@test BraketSimulator.get_variable(branched_sim, path_11, "count").val == 2

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_01 = BraketSimulator.calculate_current_state(branched_sim, path_01)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# Verify the state of qubit 2 for each path
			# Path 00: count=0, no operation, should be |0⟩
			@test abs(state_00[1]) ≈ 1.0 atol=1e-10

			# Path 01: count=1, H operation, should be (|0⟩ + |1⟩)/√2
			@test abs(state_01[3]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_01[4]) ≈ 1/sqrt(2) atol=1e-10

			# Path 10: count=1, H operation, should be (|0⟩ + |1⟩)/√2
			@test abs(state_10[5]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_10[6]) ≈ 1/sqrt(2) atol=1e-10

			# Path 11: count=2, X operation, should be |1⟩
			@test abs(state_11[8]) ≈ 1.0 atol=1e-10
		end

		@testset "4.2 Additional data types and operations" begin
			# Create an OpenQASM program with various data types
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Float data type
			float[64] rotate = 0.5;

			// Array data type
			array[int[32], 3] counts = {0, 0, 0};

			// Initialize qubits
			h q[0];
			h q[1];

			// Measure qubits
			b[0] = measure q[0];
			b[1] = measure q[1];

			// Update counts based on measurements
			if (b[0] == 1) {
				counts[0] = counts[0] + 1;
			}
			if (b[1] == 1) {
				counts[1] = counts[1] + 1;
			}
			counts[2] = counts[0] + counts[1];

			// Use float value to control rotation
			if (counts[2] > 0) {
				// Apply rotation based on angle
				U(rotate * pi, 0.0, 0.0) q[0];
			}
			"""

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 4 paths (2^2 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# Find paths for each measurement combination
			path_00 = nothing
			path_01 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 1
					path_01 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all four paths
			@test !isnothing(path_00)
			@test !isnothing(path_01)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Verify the counts array for each path
			# Path 00: counts = [0, 0, 0]
			counts_00 = BraketSimulator.get_variable(branched_sim, path_00, "counts").val
			@test counts_00[1] == 0
			@test counts_00[2] == 0
			@test counts_00[3] == 0

			# Path 01: counts = [0, 1, 1]
			counts_01 = BraketSimulator.get_variable(branched_sim, path_01, "counts").val
			@test counts_01[1] == 0
			@test counts_01[2] == 1
			@test counts_01[3] == 1

			# Path 10: counts = [1, 0, 1]
			counts_10 = BraketSimulator.get_variable(branched_sim, path_10, "counts").val
			@test counts_10[1] == 1
			@test counts_10[2] == 0
			@test counts_10[3] == 1

			# Path 11: counts = [1, 1, 2]
			counts_11 = BraketSimulator.get_variable(branched_sim, path_11, "counts").val
			@test counts_11[1] == 1
			@test counts_11[2] == 1
			@test counts_11[3] == 2

			# Verify that the rotation was applied for paths 01, 10, and 11
			# For path 00, no rotation should be applied

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_01 = BraketSimulator.calculate_current_state(branched_sim, path_01)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# For path_00, no rotation, qubit 0 should remain in the state it was measured in (|0⟩)
			@test abs(state_00[1]) ≈ 1.0 atol=1e-10

			# For paths 01, 10, and 11, a rotation was applied
			# The rotation is u3(0.5*pi, 0, 0) which is a rotation around the X-axis by 0.5*pi
			# This transforms |0⟩ to (|0⟩ + |1⟩)/√2 and |1⟩ to (|0⟩ - |1⟩)/√2
			# However, since the qubits are measured, we need to check the specific states
			# We'll skip detailed state verification for these cases
		end

		@testset "4.3 Type casting operations" begin
			# Create an OpenQASM program with type casting
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Initialize variables of different types
			int[32] int_val = 3;
			float[64] float_val = 2.5;

			// Type casting
			int[32] truncated_float = int(float_val);  // Should be 2
			float[64] float_from_int = float(int_val);  // Should be 3.0

			// Use bit casting
			bit[32] bits_from_int = bit[32](int_val);  // Binary representation of 3
			int[32] int_from_bits = int[32](bits_from_int);  // Should be 3 again

			// Initialize qubits based on casted values
			h q[0];
			h q[1];

			// Measure qubits
			b[0] = measure q[0];
			b[1] = measure q[1];

			// Use casted values in conditionals
			if (b[0] == 1 && truncated_float == 2) {
				// Apply X to qubit 0 if b[0]=1 and truncated_float=2
				x q[0];
			}

			if (b[1] == 1 && int_from_bits == 3) {
				// Apply Z to qubit 1 if b[1]=1 and int_from_bits=3
				z q[1];
			}
			"""

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 4 paths (2^2 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# Find paths for each measurement combination
			path_00 = nothing
			path_01 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 1
					path_01 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all four paths
			@test !isnothing(path_00)
			@test !isnothing(path_01)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Verify the casted values
			for path_idx in branched_sim.active_paths
				# Check truncated_float
				truncated_float = BraketSimulator.get_variable(branched_sim, path_idx, "truncated_float").val
				@test truncated_float == 2

				# Check float_from_int
				float_from_int = BraketSimulator.get_variable(branched_sim, path_idx, "float_from_int").val
				@test float_from_int ≈ 3.0 atol=1e-10

				# Check int_from_bits
				int_from_bits = BraketSimulator.get_variable(branched_sim, path_idx, "int_from_bits").val
				@test int_from_bits == 3
			end

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_01 = BraketSimulator.calculate_current_state(branched_sim, path_01)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# For path_10, b[0]=1 and truncated_float=2, so X should be applied to qubit 0
			# This means qubit 0 should be in |0⟩ state (since it was measured as |1⟩ and then flipped)
			@test abs(state_10[1]) ≈ 1.0 atol=1e-10

			# For path_11, b[0]=1 and truncated_float=2, so X should be applied to qubit 0
			# Also, b[1]=1 and int_from_bits=3, so Z should be applied to qubit 1
			# Qubit 0 should be in |0⟩ state (since it was measured as |1⟩ and then flipped)
			# Qubit 1 should be in |1⟩ state (since it was measured as |1⟩ and Z doesn't change the basis)
			@test abs(state_11[2]) ≈ 1.0 atol=1e-10

			# For path_01, b[1]=1 and int_from_bits=3, so Z should be applied to qubit 1
			# Qubit 1 should be in |1⟩ state (since it was measured as |1⟩ and Z doesn't change the basis)
			@test abs(state_01[2]) ≈ 1.0 atol=1e-10
		end
		@testset "4.4 Complex Classical Operations" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[3] q;
			bit[3] b;
			int[32] x = 5;
			float[64] y = 2.5;

			// Arithmetic operations
			float[64] w = y / 2.0;

			// Bitwise operations
			int[32] z = x * 2 + 3;
			int[32] bit_ops = (x << 1) | 3;

			h q[0];
			if (z > 10) {
				x q[1];
			}
			if (w < 2.0) {
				z q[2];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(3, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify classical variable computations
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "z").val == 13
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "w").val ≈ 1.25
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "bit_ops").val == 11
		end
	end

	@testset "5. Control Flow Structures" begin
		@testset "5.1 Loop dependent on measurement results" begin
			# Create an OpenQASM program with a loop dependent on measurement results
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

			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 4 paths:
			# 1. Measured 1 on first attempt
			# 2. Measured 0 on first attempt, then 1 on second attempt
			# 3. Measured 0 on first and second attempts, then 1 on third attempt
			# 4. Measured 0 on all three attempts
			@test length(branched_sim.active_paths) == 4

			# Find paths for each count value
			paths_by_count = Dict{Int, Vector{Int}}()

			for path_idx in branched_sim.active_paths
				count = BraketSimulator.get_variable(branched_sim, path_idx, "count").val
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
					b = BraketSimulator.get_variable(branched_sim, path_idx, "b").val
					@test b == 1

					# Calculate the final state for this path
					state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

					# Qubit 1 should be in state |1⟩
					@test abs(state[2]) ≈ 1.0 atol=1e-10 || abs(state[4]) ≈ 1.0 atol=1e-10
				end
			end

			# For paths with count = 3, check each path individually
			for path_idx in paths_by_count[3]
				b = BraketSimulator.get_variable(branched_sim, path_idx, "b").val

				# Calculate the final state for this path
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

				if b == 1
					# If b=1, qubit 1 should be in state |1⟩
					@test abs(state[2]) ≈ 1.0 atol=1e-10 || abs(state[4]) ≈ 1.0 atol=1e-10
				else
					# If b=0, qubit 1 should be in state |0⟩
					@test abs(state[1]) ≈ 1.0 atol=1e-10 || abs(state[3]) ≈ 1.0 atol=1e-10
				end
			end
		end

		@testset "5.2 For loop operations" begin
			# Create an OpenQASM program with a for loop
			qasm_source = """
			OPENQASM 3.0;
			qubit[4] q;
			bit[4] b;
			int[32] sum = 0;

			// Initialize all qubits to |+⟩ state
			for uint i in [0:3] {
				h q[i];
			}

			// Measure all qubits
			for uint i in [0:3] {
				b[i] = measure q[i];
			}

			// Count the number of 1s measured
			for uint i in [0:3] {
				if (b[i] == 1) {
					sum = sum + 1;
				}
			}

			// Apply operations based on the sum
			switch (sum) {
			case 0 {

			}
			case 1 {
				x q[0];  // Apply X to qubit 0
			}
			case 2 {
				h q[0];  // Apply H to qubit 0
			}
			case 3 {
				z q[0];  // Apply Z to qubit 0
			}
			default {
				y q[0];  // Apply Y to qubit 0
			}
			}
			"""

			# Create a simulator
			simulator = StateVectorSimulator(4, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 16 paths (2^4 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 16

			# Group paths by sum value
			paths_by_sum = Dict{Int, Vector{Int}}()

			for path_idx in branched_sim.active_paths
				sum_val = BraketSimulator.get_variable(branched_sim, path_idx, "sum").val
				if !haskey(paths_by_sum, sum_val)
					paths_by_sum[sum_val] = Int[]
				end
				push!(paths_by_sum[sum_val], path_idx)
			end

			# Verify that we have paths for all possible sums (0 to 4)
			for sum_val in 0:4
				@test haskey(paths_by_sum, sum_val)
			end

			# Verify the number of paths for each sum value
			# There should be binomial(4, k) paths with sum = k
			@test length(paths_by_sum[0]) == 1   # 1 way to get sum=0 (all 0s)
			@test length(paths_by_sum[1]) == 4   # 4 ways to get sum=1 (one 1, three 0s)
			@test length(paths_by_sum[2]) == 6   # 6 ways to get sum=2 (two 1s, two 0s)
			@test length(paths_by_sum[3]) == 4   # 4 ways to get sum=3 (three 1s, one 0)
			@test length(paths_by_sum[4]) == 1   # 1 way to get sum=4 (all 1s)

			# Check the state of qubit 0 for each sum value
			for (sum_val, path_indices) in paths_by_sum
				for path_idx in path_indices
					state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

					if sum_val == 0
						# No operation, qubit 0 should be in the state it was measured in
						# Since we're checking qubit 0 specifically, we need to check if it was measured as 0 or 1
						measured_val = branched_sim.measurements[path_idx]["q[0]"][1]
						if measured_val == 0
							# Should be in |0⟩ state
							@test abs(state[1]) ≈ 1.0 atol=1e-10
						else
							# Should be in |1⟩ state
							@test abs(state[2]) ≈ 1.0 atol=1e-10
						end
					elseif sum_val == 1
						# X operation, qubit 0 should be in |1⟩ state
						@test abs(state[2]) ≈ 1.0 atol=1e-10 || abs(state[4]) ≈ 1.0 atol=1e-10 ||
																						 abs(state[6]) ≈ 1.0 atol=1e-10 || abs(state[8]) ≈ 1.0 atol=1e-10 ||
																																					abs(state[10]) ≈ 1.0 atol=1e-10 || abs(state[12]) ≈ 1.0 atol=1e-10 ||
																																																				 abs(state[14]) ≈ 1.0 atol=1e-10 ||
																																																										   abs(state[16]) ≈ 1.0 atol=1e-10
					elseif sum_val == 2
						# H operation, qubit 0 should be in (|0⟩ + |1⟩)/√2 state
						# This is harder to test directly due to the entanglement with other qubits
						# We'll skip detailed state verification for this case
					elseif sum_val == 3
						# Z operation, doesn't change basis states, so we can't easily verify
						# We'll skip detailed state verification for this case
					elseif sum_val == 4
						# Y operation, qubit 0 should be in i|1⟩ state if it was |0⟩, or -i|0⟩ if it was |1⟩
						# This is harder to test directly due to the entanglement with other qubits
						# We'll skip detailed state verification for this case
					end
				end
			end
		end

		@testset "5.3 Complex Control Flow" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;
			int[32] count = 0;

			while (count < 2) {
				h q[count];
				b[count] = measure q[count];
				if (b[count] == 1) {
					break;
				}
				count = count + 1;
			}

			// Apply operations based on final count
			if (count == 0) {
				x q[1];
			} elif (count == 1) {
				z q[1];
			} else {
				h q[1];
			}
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify branching paths exist
			@test length(branched_sim.active_paths) > 0
		end

		@testset "5.4 Array Operations and Indexing" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[4] q;
			bit[4] b;
			array[int[32], 4] arr = {1, 2, 3, 4};

			// Array operations
			for uint i in [0:3] {
				if (arr[i] % 2 == 0) {
					h q[i];
				}
			}

			// Measure all qubits
			for uint i in [0:3] {
				b[i] = measure q[i];
			}
			"""

			simulator = StateVectorSimulator(4, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have 16 paths (2^4 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 16
		end
	end

	@testset "6. Advanced Quantum Circuits" begin
		@testset "6.1 Quantum teleportation" begin
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

			// At this point, qubit 2 should be in the |+⟩ state
			"""

			# Create a simulator
			simulator = StateVectorSimulator(3, 1000)

			# Evolve the program using the branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# We should have 4 paths (one for each combination of measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# Find paths for each measurement combination
			path_00 = nothing
			path_01 = nothing
			path_10 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 1
					path_01 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 0
					path_10 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found all four paths
			@test !isnothing(path_00)
			@test !isnothing(path_01)
			@test !isnothing(path_10)
			@test !isnothing(path_11)

			# Calculate the final states for all paths
			state_00 = BraketSimulator.calculate_current_state(branched_sim, path_00)
			state_01 = BraketSimulator.calculate_current_state(branched_sim, path_01)
			state_10 = BraketSimulator.calculate_current_state(branched_sim, path_10)
			state_11 = BraketSimulator.calculate_current_state(branched_sim, path_11)

			# For all paths, qubit 2 should be in the |+⟩ state (|0⟩ + |1⟩)/√2
			# regardless of the measurement outcomes, due to the corrections

			# Extract the reduced state of qubit 2 for each path
			# Since qubits 0 and 1 are measured, we need to check the state of qubit 2 only

			# For path_00, no corrections needed
			@test abs(state_00[1]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_00[2]) ≈ 1/sqrt(2) atol=1e-10

			# For path_01, X correction
			@test abs(state_01[3]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_01[4]) ≈ 1/sqrt(2) atol=1e-10

			# For path_10, Z correction
			@test abs(state_10[5]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_10[6]) ≈ 1/sqrt(2) atol=1e-10

			# For path_11, X and Z corrections
			@test abs(state_11[7]) ≈ 1/sqrt(2) atol=1e-10
			@test abs(state_11[8]) ≈ 1/sqrt(2) atol=1e-10
		end

		@testset "6.2 Quantum Phase Estimation" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[4] q;  // 3 counting qubits + 1 eigenstate qubit
			bit[3] b;

			// Initialize eigenstate qubit
			x q[3];

			// Apply QFT
			for uint i in [0:2] {
				h q[i];
			}

			// Controlled phase rotations
			ctrl @ phase(pi/2) q[0], q[3];
			ctrl @ phase(pi/4) q[1], q[3];
			ctrl @ phase(pi/8) q[2], q[3];

			// Inverse QFT
			for uint i in [2:-1:0] {
				for uint j in [i-1:-1:0] {
					ctrl @ phase(-pi/float(2^(i-j))) q[j], q[i];
				}
				h q[i];
			}

			// Measure counting qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			b[2] = measure q[2];
			"""

			simulator = StateVectorSimulator(4, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify expected number of paths
			@test length(branched_sim.active_paths) == 8  # 2^3 possible measurement outcomes
		end

		@testset "6.3 Dynamic Circuit Features" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;
			float[64] angle = pi/4;

			// Dynamic rotation angles
			rx(angle) q[0];
			ry(angle*2) q[1];

			b[0] = measure q[0];

			// Dynamic phase based on measurement
			if (b[0] == 1) {
				angle = angle * 2;
			} else {
				angle = angle / 2;
			}

			rz(angle) q[1];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify branching based on measurement
			@test length(branched_sim.active_paths) == 4

			# Verify different angle values in different paths
			for path_idx in branched_sim.active_paths
				angle = BraketSimulator.get_variable(branched_sim, path_idx, "angle").val
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				if b0 == 1
					@test angle ≈ pi/2 atol=1e-10
				else
					@test angle ≈ pi/8 atol=1e-10
				end
			end
		end
	end

	@testset "7. Custom Quantum Components" begin
		@testset "7.1 Custom Gates and Subroutines" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Define custom gate
			def custom_gate(qubit q) {
				h q;
				t q;
				h q;
			}

			// Define subroutine
			def measure_and_reset(qubit q, bit b) -> bit {
				b = measure q;
				if (b == 1) {
					x q;
				}
				return b;
			}

			custom_gate(q[0]);
			b[0] = measure_and_reset(q[1], b[1]);
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify custom gate application
			@test length(branched_sim.instruction_sequences[1]) > 0
		end
	end

	@testset "8. Error Handling and Edge Cases" begin
		@testset "8.1 Maximum Recursion" begin
			# Create a circuit with deeply nested conditionals
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;
			int[32] depth = 0;

			h q[0];
			b[0] = measure q[0];

			while (depth < 10) {
				if (b[0] == 1) {
					h q[1];
					b[1] = measure q[1];
					if (b[1] == 1) {
						x q[0];
						b[0] = measure q[0];
					}
				}
				depth = depth + 1;
			}
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			@test length(branched_sim.active_paths) > 0
		end
	end
end
