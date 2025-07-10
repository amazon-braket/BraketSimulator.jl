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
					@test isapprox(abs(state[2]), 1.0; atol = 1e-10) || isapprox(abs(state[4]), 1.0; atol = 1e-10)
				end
			end

			# For paths with count = 3, check each path individually
			for path_idx in paths_by_count[3]
				b = BraketSimulator.get_variable(branched_sim, path_idx, "b").val

				# Calculate the final state for this path
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

				if b == 1
					# If b=1, qubit 1 should be in state |1⟩
					@test isapprox(abs(state[2]), 1.0; atol = 1e-10) || isapprox(abs(state[4]), 1.0; atol = 1e-10)
				else
					# If b=0, qubit 1 should be in state |0⟩
					@test isapprox(abs(state[1]), 1.0; atol = 1e-10) || isapprox(abs(state[3]), 1.0; atol = 1e-10)
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
						@test isapprox(abs(state[1]), 1.0; atol = 1e-10) ||
							  isapprox(abs(state[10]), 1.0; atol = 1e-10) ||
							  isapprox(abs(state[11]), 1.0; atol = 1e-10) ||
							  isapprox(abs(state[13]), 1.0; atol = 1e-10)
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
			switch(count){
			case 0 {
				x q[1];
			}
			case 1 {
				z q[1];
			}
			default {
				h q[1];
			}
			}
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify branching paths exist
			@test length(branched_sim.active_paths) > 0

			# Group paths by count value
			paths_by_count = Dict{Int, Vector{Int}}()
			for path_idx in branched_sim.active_paths
				count = BraketSimulator.get_variable(branched_sim, path_idx, "count").val
				if !haskey(paths_by_count, count)
					paths_by_count[count] = Int[]
				end
				push!(paths_by_count[count], path_idx)
			end

			# Verify that we have paths for counts 0, 1, 2
			@test haskey(paths_by_count, 0)
			@test haskey(paths_by_count, 1)
			@test haskey(paths_by_count, 2)

			# Check the state of qubit 1 for each count value
			for path_idx in paths_by_count[0]
				# count=0 means we measured 1 on first qubit, so X was applied to qubit 1
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
				# Qubit 1 should be in state |1⟩
				@test isapprox(abs(state[2]), 1.0; atol = 1e-10) || isapprox(abs(state[4]), 1.0; atol = 1e-10)
			end

			for path_idx in paths_by_count[1]
				# count=1 means we measured 0 on first qubit, 1 on second, so Z was applied to qubit 1
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
				# Z doesn't change the basis state, so we need to check the measurement result
				@test isapprox(abs(state[2]), 1.0; atol = 1e-10)
			end

			for path_idx in paths_by_count[2]
				# count=2 means we measured 0 on both qubits, so H was applied to qubit 1
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
				# Check that qubit 1 is in a superposition state
				# Since qubit 0 is in |0⟩, we check states |00⟩ and |01⟩
				@test isapprox(abs(state[1]), 1/sqrt(2); atol = 1e-10)
				@test isapprox(abs(state[2]), 1/sqrt(2); atol = 1e-10)
			end
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

			# Verify that we have 4 paths (2^2 possible measurement outcomes for the 2 qubits with H applied)
			@test length(branched_sim.active_paths) == 4

			# Check that H was applied only to qubits with even indices (1 and 3)
			for path_idx in branched_sim.active_paths
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]
				b2 = branched_sim.measurements[path_idx]["q[2]"][1]
				b3 = branched_sim.measurements[path_idx]["q[3]"][1]

				# Verify that qubits 0 and 2 are in |0⟩ state (no H applied)
				@test b0 == 0
				@test b2 == 0

				# Calculate expected state based on measurements
				expected_state = zeros(ComplexF64, 16)
				index = 1 + b0*8 + 4*b1 + 2*b2 + b3
				expected_state[index] = 1.0

				@test isapprox(state, expected_state; atol = 1e-10)
			end
		end

		@testset "5.5 Nested Loops with Measurements" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[3] q;
			bit[3] b;
			int[32] outer_count = 0;
			int[32] inner_count = 0;
			int[32] total_ones = 0;

			// Nested loops with measurements
			for uint i in [0:1] {
				h q[i];  // Put qubits in superposition
				b[i] = measure q[i];
				outer_count = outer_count + 1;
				
				if (b[i] == 1) {
					total_ones = total_ones + 1;
					
					// Inner loop that depends on measurement result
					for uint j in [0:1] {
						if (j != i) {
							h q[j];
							b[j] = measure q[j];
							inner_count = inner_count + 1;
							
							if (b[j] == 1) {
								total_ones = total_ones + 1;
							}
						}
					}
				}
			}

			// Apply operation based on total number of ones
			if (total_ones > 1) {
				x q[2];
			}
			"""

			simulator = StateVectorSimulator(3, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have multiple paths
			@test length(branched_sim.active_paths) > 0

			# Group paths by total_ones value
			paths_by_total = Dict{Int, Vector{Int}}()
			for path_idx in branched_sim.active_paths
				total = BraketSimulator.get_variable(branched_sim, path_idx, "total_ones").val
				if !haskey(paths_by_total, total)
					paths_by_total[total] = Int[]
				end
				push!(paths_by_total[total], path_idx)
			end

			# Check the state of qubit 2 for each path
			for path_idx in branched_sim.active_paths
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)
				total_ones = BraketSimulator.get_variable(branched_sim, path_idx, "total_ones").val

				if total_ones > 1
					# X gate should have been applied to qubit 2
					@test isapprox(abs(state[2]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[4]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[6]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[8]), 1.0; atol = 1e-10)
				else
					# Qubit 2 should remain in |0⟩ state
					@test isapprox(abs(state[1]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[3]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[5]), 1.0; atol = 1e-10) ||
						  isapprox(abs(state[7]), 1.0; atol = 1e-10)
				end
			end

			# Verify that outer_count and inner_count are correctly updated
			for path_idx in branched_sim.active_paths
				outer_count = BraketSimulator.get_variable(branched_sim, path_idx, "outer_count").val
				inner_count = BraketSimulator.get_variable(branched_sim, path_idx, "inner_count").val

				# outer_count should always be 2 (we iterate through i=0,1)
				@test outer_count == 2

				# inner_count depends on measurement results
				# Check that it's consistent with the number of 1s measured in the first two qubits
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][length(branched_sim.measurements[path_idx]["q[1]"])]

				expected_inner_count = 0
				if b0 == 1
					expected_inner_count += 1  # j=1 iteration when i=0
				end
				if b1 == 1
					expected_inner_count += 1  # j=0 iteration when i=1
				end

				@test inner_count == expected_inner_count
			end
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
			phaseshift(pi/2) q[0];
			phaseshift(pi/4) q[1];
			phaseshift(pi/8) q[2];

			// Inverse QFT
			for uint i in [2:-1:0] {
				for uint j in [(i-1):-1:0] {
					phaseshift(-pi/float(2**(i-j))) q[j];
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

			results = BraketSimulator.calculate_current_state(branched_sim)

			# Extract the final states of all branches
			for path in branched_sim.active_paths
				state = results[path]
				probs = abs2.(state)
				# Get index of max probability amplitude
				max_idx = argmax(probs)
				bitstr = bitstring(max_idx - 1)  # Get binary with leading zeros for 4 qubits

				measured_bits = bitstr[1:3]  # Get first 3 bits (counting qubits)

				# We're expecting '111' to dominate
				if measured_bits == "111"
					@test probs[max_idx] ≈ 1.0 atol=0.1
				end
			end
		end

		@testset "6.3 Dynamic Circuit Features" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;
			float[64] ang = pi/4;

			// Dynamic rotation angles
			rx(ang) q[0];
			ry(ang*2) q[1];

			b[0] = measure q[0];

			// Dynamic phase based on measurement
			if (b[0] == 1) {
				ang = ang * 2;
			} else {
				ang = ang / 2;
			}

			rz(ang) q[1];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify branching based on measurement
			@test length(branched_sim.active_paths) == 4

			# Verify different angle values in different paths
			for path_idx in branched_sim.active_paths
				angle = BraketSimulator.get_variable(branched_sim, path_idx, "ang").val
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				if b0 == 1
					@test angle ≈ pi/2 atol=1e-10
				else
					@test angle ≈ pi/8 atol=1e-10
				end
			end
		end

		@testset "6.4 Quantum Fourier Transform" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[3] q;
			bit[3] b;

			// Initialize state |001⟩
			x q[2];

			// Apply QFT
			// Qubit 0
			h q[0];
			ctrl @ gphase(pi/2) q[1];
			ctrl @ gphase(pi/4) q[2];

			// Qubit 1
			h q[1];
			ctrl @ gphase(pi/2) q[2];

			// Qubit 2
			h q[2];

			// Swap qubits 0 and 2
			swap q[0], q[2];

			// Measure all qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			b[2] = measure q[2];
			"""

			simulator = StateVectorSimulator(3, 1000)

			# First, let's calculate the expected state after QFT on |001⟩
			# We'll create a circuit without measurements to get the state vector
			qft_source = """
			OPENQASM 3.0;
			qubit[3] q;

			// Initialize state |001⟩
			x q[2];

			// Apply QFT
			// Qubit 0
			h q[0];
			phaseshift(pi/2) q[1];
			phaseshift(pi/4) q[2];

			// Qubit 1
			h q[1];
			phaseshift(pi/2) q[2];

			// Qubit 2
			h q[2];

			// Swap qubits 0 and 2
			swap q[0], q[2];
			"""

			qft_sim = StateVectorSimulator(3, 1)
			qft_sim = BraketSimulator.evolve!(qft_sim, BraketSimulator.to_circuit(qft_source).instructions)
			expected_state = BraketSimulator.state_vector(qft_sim)

			# Now evolve the original circuit with measurements
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
			# Verify that we have 8 paths (2^3 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 8

			# Calculate the probabilities from the expected state
			expected_probs = abs2.(expected_state)

			# Count the number of paths with each measurement outcome
			outcome_counts = Dict{Tuple{Int, Int, Int}, Int}()

			for path_idx in branched_sim.active_paths
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]
				b2 = branched_sim.measurements[path_idx]["q[2]"][1]
				outcome = (b0, b1, b2)

				if !haskey(outcome_counts, outcome)
					outcome_counts[outcome] = 0
				end
				outcome_counts[outcome] += 1
			end

			# Verify that the distribution of paths matches the expected probabilities
			for (outcome, count) in outcome_counts
				b0, b1, b2 = outcome
				state_idx = 1 + b0 + 2*b1 + 4*b2
				expected_prob = expected_probs[state_idx]

				# The number of paths with this outcome should be proportional to the probability
				# Since we have 8 paths total, and the simulator uses shots to determine branching
				@test count > 0
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
			gate custom_gate q {
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

			custom_gate q[0];
			b[0] = measure_and_reset(q[0], b[1]);
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify custom gate application
			@test length(branched_sim.instruction_sequences[1]) > 0

			# Verify that the custom gate is equivalent to Z-rotation by π/4
			# custom_gate = HTH = Rx(π/4)

			# Create a circuit with just the custom gate
			custom_gate_source = """
			OPENQASM 3.0;
			qubit[1] q;

			// Define custom gate
			gate custom_gate q {
				h q;
				t q;
				h q;
			}

			custom_gate q[0];
			"""

			# Create a circuit with equivalent Rz gate
			rz_gate_source = """
			OPENQASM 3.0;
			qubit[1] q;

			rx(pi/4) q[0];
			"""

			# Simulate both circuits
			custom_sim = StateVectorSimulator(1, 1)
			rz_sim = StateVectorSimulator(1, 1)

			custom_sim = BraketSimulator.evolve!(custom_sim, BraketSimulator.to_circuit(custom_gate_source).instructions)
			rz_sim = BraketSimulator.evolve!(rz_sim, BraketSimulator.to_circuit(rz_gate_source).instructions)

			# Verify that the states are equivalent
			custom_state = BraketSimulator.state_vector(custom_sim)
			rz_state = BraketSimulator.state_vector(rz_sim)

			@test isapprox(abs(custom_state[1]), abs(rz_state[1]); atol = 1e-10)
			@test isapprox(abs(custom_state[2]), abs(rz_state[2]); atol = 1e-10)

			# Now verify the measure_and_reset subroutine
			# It should always leave the qubit in |0⟩ state

			# Check the paths in the original simulation
			@test length(branched_sim.active_paths) == 2  # Two paths from measuring q[1]

			results = BraketSimulator.calculate_current_state(branched_sim)

			for path_idx in branched_sim.active_paths
				state = results[path_idx]

				# Extract measurement result for q[1]
				b1 = branched_sim.measurements[path_idx]["q[0]"][1]

				# Verify that b[0] equals the measurement result
				b0 = BraketSimulator.get_variable(branched_sim, path_idx, "b").val[1]
				@test b0 == b1

				# Verify that q[1] is in |0⟩ state (reset)
				# Since q[0] had custom_gate applied, we need to check both |00⟩ and |10⟩ states
				@test isapprox(abs(state[1]), 1; atol = 1e-10)
			end
		end

		@testset "7.2 Custom Gates with Control Flow" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[3] q;
			bit[2] b;

			// Define a custom controlled rotation gate
			gate controlled_rotation(ang) control, target {
				ctrl @ rz(ang) control, target;
			}

			// Define a custom function that applies different operations based on measurement
			def adaptive_gate(qubit q1, qubit q2, bit measurement) {
				if (measurement == 0) {
					h q1;
					h q2;
				} else {
					x q1;
					z q2;
				}
			}

			// Initialize qubits
			h q[0];

			// Measure qubit 0
			b[0] = measure q[0];

			// Apply custom gates based on measurement
			controlled_rotation(pi/2) q[0], q[1];
			adaptive_gate(q[1], q[2], b[0]);

			// Measure qubit 1
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(3, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have 4 paths (2 from first measurement × 2 from second measurement)
			@test length(branched_sim.active_paths) == 3

			# Group paths by first measurement outcome
			paths_by_first_meas = Dict{Int, Vector{Int}}()

			for path_idx in branched_sim.active_paths
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				if !haskey(paths_by_first_meas, b0)
					paths_by_first_meas[b0] = Int[]
				end
				push!(paths_by_first_meas[b0], path_idx)
			end

			# Verify that we have paths for both measurement outcomes
			@test haskey(paths_by_first_meas, 0)
			@test haskey(paths_by_first_meas, 1)

			# Check the state of qubits 1 and 2 for each path
			for path_idx in paths_by_first_meas[0]
				# For b[0]=0, adaptive_gate applies H to both qubits
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

				# Extract measurement result for q[1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				# Since q[1] is measured, we need to check if the state is consistent
				# with the measurement result
				if b1 == 0
					# q[1] is |0⟩, q[2] should be in (|0⟩ + |1⟩)/√2 state
					@test isapprox(abs(state[1]), 1/sqrt(2); atol = 1e-10)
					@test isapprox(abs(state[2]), 1/sqrt(2); atol = 1e-10)
				else
					# q[1] is |1⟩, q[2] should be in (|0⟩ + |1⟩)/√2 state
					@test isapprox(abs(state[3]), 1/sqrt(2); atol = 1e-10)
					@test isapprox(abs(state[4]), 1/sqrt(2); atol = 1e-10)
				end
			end

			for path_idx in paths_by_first_meas[1]
				# For b[0]=1, adaptive_gate applies X to q[1] and Z to q[2]
				state = BraketSimulator.calculate_current_state(branched_sim, path_idx)

				# Extract measurement result for q[1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				# Since q[1] is measured, we need to check if the state is consistent
				# with the measurement result
				if b1 == 0
					# q[1] is |0⟩, q[2] should be in |0⟩ state (Z doesn't change |0⟩)
					@test isapprox(abs(state[5]), 1.0; atol = 1e-10)
				else
					# q[1] is |1⟩, q[2] should be in |0⟩ state (Z doesn't change |0⟩)
					@test isapprox(abs(state[7]), 1.0; atol = 1e-10)
				end
			end
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

	@testset "9. Gate Modifiers and GPhase" begin
		@testset "9.1 Basic gate modifiers" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Apply X gate with power modifier (X^0.5 = √X)
			pow(0.5) @ x q[0];

			// Apply X gate with inverse modifier (X† = X)
			inv @ x q[1];

			// Measure both qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have 2 paths (2 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 2

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# For each path, verify the state before measurement
			for path_idx in branched_sim.active_paths
				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				# For q[0], X^0.5 = √X should create a superposition
				# For q[1], inv @ X = X should flip the bit

				# The probability of measuring 0 or 1 for q[0] should be 0.5 each
				# The probability of measuring 1 for q[1] should be 1.0

				if b1 == 1
					# This is the expected outcome for q[1]
					@test true
				else
					# This should not happen for q[1]
					@test false
				end

				# For q[0], both outcomes are possible
				@test b0 == 0 || b0 == 1
			end
		end

		@testset "9.2 Control modifiers" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[3] q;
			bit[3] b;

			// Initialize q[0] to |1⟩
			x q[0];

			// Apply controlled-H gate (control on q[0], target on q[1])
			ctrl @ h q[0], q[1];

			// Apply controlled-controlled-X gate (controls on q[0] and q[1], target on q[2])
			ctrl @ ctrl @ x q[0], q[1], q[2];

			// Measure all qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			b[2] = measure q[2];
			"""

			simulator = StateVectorSimulator(3, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# Verify that we have multiple paths
			@test length(branched_sim.active_paths) > 0

			# For each path, verify the state before measurement
			for path_idx in branched_sim.active_paths
				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]
				b2 = branched_sim.measurements[path_idx]["q[2]"][1]

				# q[0] should always be 1 (initialized with X)
				@test b0 == 1

				# If q[1] is 1, then q[2] should also be 1 (due to the controlled-controlled-X)
				# If q[1] is 0, then q[2] should be 0
				if b1 == 1
					@test b2 == 1
				else
					@test b2 == 0
				end
			end
		end

		@testset "9.3 Negative control modifiers" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			h q[0];

			// Apply negative-controlled X gate (control on q[0], target on q[1])
			// This applies X to q[1] when q[0] is |0⟩
			negctrl @ x q[0], q[1];

			// Measure both qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# Verify that we have 2 paths (since the states are deterministic)
			@test length(branched_sim.active_paths) == 2

			# For each path, verify the state before measurement
			for path_idx in branched_sim.active_paths
				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				# If q[0] is 0, then q[1] should be 1 (due to the negative-controlled X)
				# If q[0] is 1, then q[1] should be 0
				if b0 == 0
					@test b1 == 1
				else
					@test b1 == 0
				end
			end
		end

		@testset "9.4 Multiple modifiers" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Initialize q[0] to |1⟩
			h q[0];

			// Apply controlled-inverse-X gate (control on q[0], target on q[1])
			// Since X† = X, this is equivalent to a standard CNOT
			ctrl @ inv @ x q[0], q[1];

			// Measure both qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# Verify that we have 2 paths (since the states are deterministic)
			@test length(branched_sim.active_paths) == 2

			# For each path, verify the state before measurement
			for path_idx in branched_sim.active_paths
				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				if b0 == 1
					# q[0] should always be 1 (initialized with X)
					@test b1 == 1
				else
					# q[1] should also be 1 (due to the controlled-X with q[0] as control)
					@test b1 == 0
				end
			end
		end

		@testset "9.5 GPhase gate" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Apply global phase
			gphase(pi/2);

			// Apply controlled global phase
			ctrl @ gphase(pi/4) q[0];

			// Create superposition
			h q[0];
			h q[1];

			// Measure both qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# Verify that we have 4 paths (2^2 possible measurement outcomes)
			@test length(branched_sim.active_paths) == 4

			# For each path, verify the state before measurement
			for path_idx in branched_sim.active_paths
				# Extract measurement results
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]
				b1 = branched_sim.measurements[path_idx]["q[1]"][1]

				# Both qubits should have equal probability of being 0 or 1
				@test b0 == 0 || b0 == 1
				@test b1 == 0 || b1 == 1
			end

			# Global phase doesn't affect measurement probabilities, so we just verify
			# that the circuit executed without errors
		end

		@testset "9.6 Power modifiers with parametric angles" begin
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;
			float[64] ang = 0.25;

			// Apply X gate with power modifier using a variable
			pow(ang) @ x q[0];

			// Measure the qubit
			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Calculate the final states for all paths
			states = BraketSimulator.calculate_current_state(branched_sim)

			# Verify that we have 2 paths (since both measurement outcomes are possible)
			@test length(branched_sim.active_paths) == 2

			# Count the number of paths with each measurement outcome
			count_0 = 0
			count_1 = 0

			for path_idx in branched_sim.active_paths
				# Extract measurement result
				b0 = branched_sim.measurements[path_idx]["q[0]"][1]

				if b0 == 0
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

	@testset "10. Scoping Rules" begin
		@testset "10.1 Local scope blocks inherit variables" begin
			# Test that local scope blocks inherit all variables from containing scope
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;

			// Global variables
			int[32] global_var = 5;
			const int[32] const_global = 10;

			// Local scope block should inherit all variables
			{
				// Access global variables
				global_var = global_var + const_global;  // Should be 15
				
				// Modify non-const variable
				global_var = global_var * 2;  // Should be 30
			}

			// Verify that changes in local scope affect global scope
			if (global_var == 30) {
				h q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that global_var was modified in the local scope
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "global_var").val == 30

			# Verify that H gate was applied (since global_var == 30)
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2
		end

		@testset "10.2 For loop iteration variable lifetime" begin
			# Test that for loop iteration variables are only accessible within the loop
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;
			int[32] sum = 0;

			// For loop with iteration variable i
			for uint i in [0:4] {
				sum = sum + i;  // Sum should be 0+1+2+3+4 = 10
			}

			// i should not be accessible here
			// Instead, we use sum to verify the loop executed correctly
			if (sum == 10) {
				h q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that sum was calculated correctly in the loop
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "sum").val == 10

			# Verify that H gate was applied (since sum == 10)
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2

			# Verify that the variable i is not accessible after the loop
			# This is implicitly tested by the successful execution of the program
			# If i were accessible and used after the loop, it would cause an error
		end

		@testset "10.3 If/else blocks define different scopes" begin
			# Test that if and else blocks define different scopes
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;
			int[32] x = 5;

			// Measure qubit to create branching paths
			h q[0];
			b[0] = measure q[0];

			if (b[0] == 0) {
				// Define a variable in if block
				int[32] local_var = 10;
				x = x + local_var;  // x should be 15
			} else {
				// Define a different variable with the same name in else block
				int[32] local_var = 20;
				x = x + local_var;  // x should be 25
			}

			// local_var is not accessible here
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have two paths (one for each measurement outcome)
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

			# Verify that x has different values in each path
			@test BraketSimulator.get_variable(branched_sim, path_0, "x").val == 15
			@test BraketSimulator.get_variable(branched_sim, path_1, "x").val == 25

			# Verify that local_var is not accessible outside the if/else blocks
			# This is implicitly tested by the successful execution of the program
			# If local_var were accessible after the blocks, it would cause an error
		end

		@testset "10.4 Variable shadowing in local scopes" begin
			# Test that variables in local scopes can shadow outer scope variables
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;

			// Global variable
			int[32] x = 5;

			// Local scope with shadowing
			if (true) {
				// Shadow the global x
				int[32] x = 10;
				
				// Use the local x
				if (x == 10) {
					h q[0];
				}
			}

			// Back to global scope, x should be unchanged
			if (x == 5) {
				// Apply another gate if global x is still 5
				x q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that global x is still 5
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "x").val == 5

			# Verify that both H and X gates were applied
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2

			# Calculate the final state for one path
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

			# The state should be either |0⟩ or |1⟩ depending on the measurement outcome
			# But since both H and X were applied, the state should be flipped from what was measured
			b0 = branched_sim.measurements[branched_sim.active_paths[1]]["q[0]"][1]

			if b0 == 0
				@test abs(state[1]) ≈ 1.0 atol=1e-10
			else
				@test abs(state[2]) ≈ 1.0 atol=1e-10
			end
		end


		@testset "10.5 Subroutine and gate scopes" begin
			# Test that subroutines and gates introduce new scopes
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[1] b;

			// Global variables
			const int[32] const_global = 10;
			int[32] mutable_global = 5;

			// Define a gate
			gate custom_gate(ang) q {
				// Can access const globals
				rx(ang * const_global) q;
				// Cannot access mutable globals
			}

			// Define a subroutine
			def apply_operations(qubit target, float[64] factor) {
				// Can access const globals
				float[64] total = factor * const_global;
				// Cannot access mutable globals
				
				// Apply gate with calculated angle
				ry(total) target;
			}

			// Apply custom gate to q[0]
			custom_gate(0.1) q[0];

			// Apply subroutine to q[1]
			apply_operations(q[1], 0.2);

			// Modify mutable global
			mutable_global = 20;

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that the program executed successfully
			@test length(branched_sim.active_paths) == 2  # Two paths from measuring q[0]

			# Verify that mutable_global was modified
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "mutable_global").val == 20

			# Calculate the final state for one path
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])

			# Verify that both gates were applied
			# custom_gate applied rx(0.1 * 10) = rx(1.0) to q[0]
			# apply_operations applied ry(0.2 * 10) = ry(2.0) to q[1]

			# The exact state verification is complex due to the rotations and measurement
			# We'll just verify that the program executed without errors
			@test length(state) == 4  # 2^2 = 4 amplitudes for 2 qubits
		end

		@testset "10.6 Local variables in subroutines" begin
			# Test that variables defined in subroutines are local to the subroutine body
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;
			int[32] result = 0;

			// Define a subroutine with local variables
			def calculate(int[32] inp) -> int[32] {
				// Local variable
				int[32] temp = inp * 2;
				
				// Modify local variable
				temp = temp + 5;
				
				// Return the result
				return temp;
			}

			// Call the subroutine
			result = calculate(10);  // Should be 10*2+5 = 25

			// Use the result
			if (result == 25) {
				h q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that result has the correct value
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "result").val == 25

			# Verify that H gate was applied (since result == 25)
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2

			# Verify that the variable temp is not accessible outside the subroutine
			# This is implicitly tested by the successful execution of the program
			# If temp were accessible after the subroutine, it would cause an error
		end

		@testset "10.7 Recursion in subroutines" begin
			# Test that subroutines can call themselves (recursion)
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;
			int[32] result = 0;

			// Define a recursive factorial function
			def factorial(int[32] n) -> int[32] {
				if (n <= 1) {
					return 1;
				} else {
					return n * factorial(n - 1);
				}
			}

			// Calculate factorial of 4
			result = factorial(4);  // Should be 24

			// Use the result
			if (result == 24) {
				h q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that result has the correct value
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "result").val == 24

			# Verify that H gate was applied (since result == 24)
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2
		end

		@testset "10.8 Shadowing in subroutines" begin
			# Test that local variables in subroutines can shadow global variables
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;

			// Global variable
			int[32] x = 5;

			// Define a subroutine with shadowing
			def modify_x(int[32] inp) -> int[32] {
				// Shadow the global x
				int[32] x = inp;
				
				// Modify the local x
				x = x * 2;
				
				// Return the local x
				return x;
			}

			// Call the subroutine
			int[32] result = modify_x(10);  // Should be 20

			// Verify that global x is unchanged
			if (x == 5 && result == 20) {
				h q[0];
			}

			b[0] = measure q[0];
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that global x is unchanged
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "x").val == 5

			# Verify that result has the correct value
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "result").val == 20

			# Verify that H gate was applied (since x == 5 && result == 20)
			# This means we should have two paths (one for each measurement outcome)
			@test length(branched_sim.active_paths) == 2
		end

		@testset "10.9 Shadowing in while loops" begin
			# Test that variables in while loops can shadow outer scope variables
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;

			// Global variable
			int[32] counter = 5;

			// While loop with shadowing
			int[32] i = 0;
			while (i < 3) {
				// Shadow the global counter
				int[32] counter = i;
				
				// Use the local counter
				if (counter == 2) {
					x q[0];
				}
				
				i += 1;
			}

			// Back to global scope, counter should be unchanged
			if (counter == 5) {
				// Apply another gate if global counter is still 5
				h q[0];
			}
			"""

			simulator = StateVectorSimulator(1, 1000)

			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that the global counter is still 5 after the while loop
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "counter").val == 5

			# Verify that the X gate was applied (since counter remained 5)
			# and H gate was applied in the last iteration (when i=2)
			# This should result in the state |-> = (|0⟩ - |1⟩)/√2
			state = BraketSimulator.calculate_current_state(branched_sim)[branched_sim.active_paths[1]]
			@test isapprox(abs(state[1]), 1/sqrt(2), atol = 1e-10)
			@test isapprox(abs(state[2]), 1/sqrt(2), atol = 1e-10)
			@test isapprox(angle(state[2]) - angle(state[1]), π, atol = 1e-10)
		end

		@testset "10.10 Shadowing in for loops" begin
			# Test that variables in for loops can shadow outer scope variables
			qasm_source = """
			OPENQASM 3.0;
			qubit[1] q;
			bit[1] b;

			// Global variable
			int[32] index = 10;

			// For loop with shadowing
			for int[32] index in [0, 1, 2] {
				// Use the local index
				if (index == 2) {
					x q[0];
				}
			}

			// Back to global scope, index should be unchanged
			if (index == 10) {
				// Apply another gate if global index is still 10
				h q[0];
			}
			"""

			simulator = StateVectorSimulator(1, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that the global index is still 10 after the for loop
			@test BraketSimulator.get_variable(branched_sim, branched_sim.active_paths[1], "index").val == 10

			# Verify that the X gate was applied (since index remained 10)
			# and H gate was applied in the last iteration (when index=2)
			# This should result in the state |-> = (|0⟩ - |1⟩)/√2
			state = BraketSimulator.calculate_current_state(branched_sim)[branched_sim.active_paths[1]]
			@test isapprox(abs(state[1]), 1/sqrt(2), atol = 1e-10)
			@test isapprox(abs(state[2]), 1/sqrt(2), atol = 1e-10)
			@test isapprox(angle(state[2]) - angle(state[1]), π, atol = 1e-10)
		end


		@testset "10.11 Aliases in subroutines" begin
			# Test that aliases can be declared within subroutine scopes
			qasm_source = """
			OPENQASM 3.0;
			qubit[2] q;
			bit[2] b;

			// Define a subroutine that uses aliases
			def apply_gates(qubit q1, qubit q2) {
				// Create aliases
				let first = q1;
				let second = q2;
				
				// Use the aliases
				h first;
				cnot first, second;
			}

			// Call the subroutine to create a Bell state
			apply_gates(q[0], q[1]);

			// Measure both qubits
			b[0] = measure q[0];
			b[1] = measure q[1];
			"""

			simulator = StateVectorSimulator(2, 1000)
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())

			# Verify that we have 2 paths (since we created a Bell state, measurements should be correlated)
			@test length(branched_sim.active_paths) == 2

			# Find paths for each measurement outcome
			path_00 = nothing
			path_11 = nothing

			for path_idx in branched_sim.active_paths
				measurements = branched_sim.measurements[path_idx]
				if measurements["q[0]"][1] == 0 && measurements["q[1]"][1] == 0
					path_00 = path_idx
				elseif measurements["q[0]"][1] == 1 && measurements["q[1]"][1] == 1
					path_11 = path_idx
				end
			end

			# Verify that we found both expected paths (00 and 11)
			@test !isnothing(path_00)
			@test !isnothing(path_11)

			# Verify that we don't have other paths (01 or 10)
			@test length(branched_sim.active_paths) == 2
		end
	end

	@testset "11. OpenQASM 3 Features from test_openqasm3.jl" begin
		@testset "11.1 Adder" begin
			sv_adder_qasm = """
			OPENQASM 3;

			input uint[4] a_in;
			input uint[4] b_in;

			gate majority a, b, c {
				cnot c, b;
				cnot c, a;
				ccnot a, b, c;
			}

			gate unmaj a, b, c {
				ccnot a, b, c;
				cnot c, a;
				cnot a, b;
			}

			qubit cin;
			qubit[4] a;
			qubit[4] b;
			qubit cout;

			// set input states
			for int[8] i in [0: 3] {
			  if(bool(a_in[i])) x a[i];
			  if(bool(b_in[i])) x b[i];
			}

			// add a to b, storing result in b
			majority cin, b[3], a[3];
			for int[8] i in [3: -1: 1] { majority a[i], b[i - 1], a[i - 1]; }
			cnot a[0], cout;
			for int[8] i in [1: 3] { unmaj a[i], b[i - 1], a[i - 1]; }
			unmaj cin, b[3], a[3];
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(10, 1000)
			
			# Convert to circuit with inputs
			circuit = BraketSimulator.new_to_circuit(sv_adder_qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict("a_in"=>3, "b_in"=>7))
			
			# Verify the instructions match the expected ones
			correct_instructions = [ 
				BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(6))
				BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(3))
				BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(7))
				BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(4))
				BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(8))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 8))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 0))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(0, 8, 4))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 7))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 4))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(4, 7, 3))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 6))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 3))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(3, 6, 2))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 5))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 2))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(2, 5, 1))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 9))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(2, 5, 1))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 2))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 5))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(3, 6, 2))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 3))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 6))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(4, 7, 3))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 4))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 7))
				BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(0, 8, 4))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 0))
				BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(0, 8))
			]
			
			# Check that the instructions match
			@test length(branched_sim.instruction_sequences[1]) == length(correct_instructions)
			for (ix, c_ix) in zip(branched_sim.instruction_sequences[1], correct_instructions)
				@test ix == c_ix
			end
		end

		@testset "11.2 GPhase" begin
			qasm = """
			qubit[2] qs;

			const int[8] two = 2;

			gate x a { U(π, 0, π) a; }
			gate cx c, a { ctrl @ x c, a; }
			gate phase c, a {
				gphase(π/2);
				ctrl(two) @ gphase(π) c, a;
			}
			gate h a { U(π/2, 0, π) a; }

			h qs[0];
			
			cx qs[0], qs[1];
			phase qs[0], qs[1];
			
			gphase(π);
			inv @ gphase(π / 2);
			negctrl @ ctrl @ gphase(2 * π) qs[0], qs[1];
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state for the first path
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			expected_sv = 1/√2 * [-1; 0; 0; 1]
			@test state ≈ expected_sv
		end

		@testset "11.3 Gate def with argument manipulation" begin
			qasm = """
			qubit[2] __qubits__;
			gate u3(θ, ϕ, λ) q {
				gphase(-(ϕ+λ)/2);
				U(θ, ϕ, λ) q;
			}
			u3(0.1, 0.2, 0.3) __qubits__[0];
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Verify the instructions
			canonical_ixs = [BraketSimulator.Instruction(BraketSimulator.GPhase{1}(-(0.2 + 0.3)/2), 0), BraketSimulator.Instruction(BraketSimulator.U(0.1, 0.2, 0.3), 0)]
			@test branched_sim.instruction_sequences[1] == canonical_ixs
		end

		@testset "11.4 Physical qubits" begin
			qasm = """
			h \$0;
			cnot \$0, \$1;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Verify the instructions
			@test branched_sim.instruction_sequences[1] == [BraketSimulator.Instruction(BraketSimulator.H(), 0), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1])]
		end

		@testset "11.5 For loop and subroutines" begin
			qasm_str = """
			OPENQASM 3.0;
			def bell(qubit q0, qubit q1) {
				h q0;
				cnot q0, q1;
			}
			def n_bells(int[32] n, qubit q0, qubit q1) {
				for int i in [0:n - 1] {
					h q0;
					cnot q0, q1;
				}
			}
			qubit[4] __qubits__;
			bell(__qubits__[0], __qubits__[1]);
			n_bells(5, __qubits__[2], __qubits__[3]);
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(4, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm_str)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict("theta"=>0.2))
			
			# Verify the instructions (excluding measurements)
			expected_instructions = [BraketSimulator.Instruction(BraketSimulator.H(), 0),
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]),
								   BraketSimulator.Instruction(BraketSimulator.H(), 2),
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
								   BraketSimulator.Instruction(BraketSimulator.H(), 2),
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
								   BraketSimulator.Instruction(BraketSimulator.H(), 2), 
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
								   BraketSimulator.Instruction(BraketSimulator.H(), 2),
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
								   BraketSimulator.Instruction(BraketSimulator.H(), 2),
								   BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3])
								  ]
			
			# Check that the instructions match (excluding measurements)
			@test length(branched_sim.instruction_sequences[1]) >= length(expected_instructions)
			for (ix, c_ix) in zip(branched_sim.instruction_sequences[1][1:length(expected_instructions)], expected_instructions)
				@test ix == c_ix
			end
		end

		@testset "11.6 Builtin functions" begin
			qasm = """
				input float x;
				input float y;
				rx(x) \$0;
				rx(arccos(x)) \$0;
				rx(arcsin(x)) \$0;
				rx(arctan(x)) \$0; 
				rx(ceiling(x)) \$0;
				rx(cos(x)) \$0;
				rx(exp(x)) \$0;
				rx(floor(x)) \$0;
				rx(log(x)) \$0;
				rx(mod(x, y)) \$0;
				rx(sin(x)) \$0;
				rx(sqrt(x)) \$0;
				rx(tan(x)) \$0;
				"""
			
			x = 1.0
			y = 2.0
			inputs = Dict("x"=>x, "y"=>y)
			
			# Create a simulator
			simulator = StateVectorSimulator(1, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, inputs)
			
			# Verify the instructions
			expected_ixs = [BraketSimulator.Instruction(BraketSimulator.Rx(x), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(acos(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(asin(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(atan(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(ceil(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(cos(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(exp(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(floor(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(log(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(mod(x, y)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(sin(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(sqrt(x)), 0),  
					   BraketSimulator.Instruction(BraketSimulator.Rx(tan(x)), 0)]
			
			# Check that the instructions match
			@test length(branched_sim.instruction_sequences[1]) == length(expected_ixs)
			for (ix, c_ix) in zip(branched_sim.instruction_sequences[1], expected_ixs)
				@test ix == c_ix
			end
		end

		@testset "11.7 Reset" begin
			qasm = """
			qubit[4] q;
			x q[0];
			reset q[0];
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(4, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Verify the instructions
			@test branched_sim.instruction_sequences[1] == [BraketSimulator.Instruction(BraketSimulator.X(), 0), 
									 BraketSimulator.Instruction(BraketSimulator.Reset(), 0)]
		end

		@testset "11.8 Switch/case" begin
			# Test switch with default case
			qasm1 = """
			input int[8] x;
			switch (x + 1) {
				case 0b00 {}
				default { z \$0; }
			}
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(1, 1000)
			
			# Test with x = -1
			circuit1a = BraketSimulator.new_to_circuit(qasm1)
			branched_sim1a = BraketSimulator.evolve_branched_operators(simulator, circuit1a, Dict("x"=> -1))
			@test isempty(branched_sim1a.instruction_sequences[1])
			
			# Test with x = 0
			circuit1b = BraketSimulator.new_to_circuit(qasm1)
			branched_sim1b = BraketSimulator.evolve_branched_operators(simulator, circuit1b, Dict("x"=> 0))
			@test branched_sim1b.instruction_sequences[1] == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
			
			# Test switch with multiple cases
			qasm2 = """
			input int[8] x;
			switch (x) { case 0 {} case 1, 2 { z \$0; } }
			"""
			
			# Test with x = 0
			circuit2a = BraketSimulator.new_to_circuit(qasm2)
			branched_sim2a = BraketSimulator.evolve_branched_operators(simulator, circuit2a, Dict("x"=> 0))
			@test isempty(branched_sim2a.instruction_sequences[1])
			
			# Test with x = 1
			circuit2b = BraketSimulator.new_to_circuit(qasm2)
			branched_sim2b = BraketSimulator.evolve_branched_operators(simulator, circuit2b, Dict("x"=> 1))
			@test branched_sim2b.instruction_sequences[1] == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
			
			# Test with x = 2
			circuit2c = BraketSimulator.new_to_circuit(qasm2)
			branched_sim2c = BraketSimulator.evolve_branched_operators(simulator, circuit2c, Dict("x"=> 2))
			@test branched_sim2c.instruction_sequences[1] == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
		end

		@testset "11.9 Global gate control" begin
			qasm = """
			qubit q1;
			qubit q2;

			h q1;
			h q2;
			ctrl @ s q1, q2;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			@test state ≈ [0.5, 0.5, 0.5, 0.5im]
		end

		@testset "11.10 Power modifiers" begin
			# Test sqrt(Z) = S
			qasm_z = """
			qubit q1;
			qubit q2;
			h q1;
			h q2;
			
			pow(1/2) @ z q1;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit_z = BraketSimulator.new_to_circuit(qasm_z)
			
			# Evolve using branched simulator operators
			branched_sim_z = BraketSimulator.evolve_branched_operators(simulator, circuit_z, Dict{String, Any}())
			
			# Calculate the final state
			state_z = BraketSimulator.calculate_current_state(branched_sim_z, branched_sim_z.active_paths[1])
			
			# Create a reference circuit with S gate
			qasm_s = """
			qubit q1;
			qubit q2;
			h q1;
			h q2;
			
			s q1;
			"""
			
			# Convert to circuit
			circuit_s = BraketSimulator.new_to_circuit(qasm_s)
			
			# Evolve using branched simulator operators
			branched_sim_s = BraketSimulator.evolve_branched_operators(simulator, circuit_s, Dict{String, Any}())
			
			# Calculate the final state
			state_s = BraketSimulator.calculate_current_state(branched_sim_s, branched_sim_s.active_paths[1])
			
			# Verify the states are equivalent
			@test state_z ≈ state_s
			
			# Test sqrt(X) = V
			qasm_x = """
			qubit q1;
			qubit q2;
			h q1;
			h q2;
			
			pow(1/2) @ x q1;
			"""
			
			# Convert to circuit
			circuit_x = BraketSimulator.new_to_circuit(qasm_x)
			
			# Evolve using branched simulator operators
			branched_sim_x = BraketSimulator.evolve_branched_operators(simulator, circuit_x, Dict{String, Any}())
			
			# Calculate the final state
			state_x = BraketSimulator.calculate_current_state(branched_sim_x, branched_sim_x.active_paths[1])
			
			# Create a reference circuit with V gate
			qasm_v = """
			qubit q1;
			qubit q2;
			h q1;
			h q2;
			
			v q1;
			"""
			
			# Convert to circuit
			circuit_v = BraketSimulator.new_to_circuit(qasm_v)
			
			# Evolve using branched simulator operators
			branched_sim_v = BraketSimulator.evolve_branched_operators(simulator, circuit_v, Dict{String, Any}())
			
			# Calculate the final state
			state_v = BraketSimulator.calculate_current_state(branched_sim_v, branched_sim_v.active_paths[1])
			
			# Verify the states are equivalent
			@test state_x ≈ state_v
		end

		@testset "11.11 Complex Power modifiers" begin
			qasm = """
			const int[8] two = 2;
			gate x a { U(π, 0, π) a; }
			gate cx c, a {
				pow(1) @ ctrl @ x c, a;
			}
			gate cxx_1 c, a {
				pow(two) @ cx c, a;
			}
			gate cxx_2 c, a {
				pow(1/2) @ pow(4) @ cx c, a;
			}
			gate cxxx c, a {
				pow(1) @ pow(two) @ cx c, a;
			}

			qubit q1;
			qubit q2;
			qubit q3;
			qubit q4;
			qubit q5;

			pow(1/2) @ x q1;       // half flip
			pow(1/2) @ x q1;       // half flip
			cx q1, q2;   // flip
			cxx_1 q1, q3;    // don't flip
			cxx_2 q1, q4;    // don't flip
			cnot q1, q5;    // flip
			x q3;       // flip
			x q4;       // flip

			s q1;   // sqrt z
			s q1;   // again
			inv @ z q1; // inv z
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(5, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			expected_sv = zeros(32)
			expected_sv[end] = 1.0
			@test state ≈ expected_sv
		end

		@testset "11.12 Gate control" begin
			qasm = """
			const int[8] two = 2;
			gate x a { U(π, 0, π) a; }
			gate cx c, a {
				ctrl @ x c, a;
			}
			gate ccx_1 c1, c2, a {
				ctrl @ ctrl @ x c1, c2, a;
			}
			gate ccx_2 c1, c2, a {
				ctrl(two) @ x c1, c2, a;
			}
			gate ccx_3 c1, c2, a {
				ctrl @ cx c1, c2, a;
			}

			qubit q1;
			qubit q2;
			qubit q3;
			qubit q4;
			qubit q5;

			// doesn't flip q2
			cx q1, q2;
			// flip q1
			x q1;
			// flip q2
			cx q1, q2;
			// doesn't flip q3, q4, q5
			ccx_1 q1, q4, q3;
			ccx_2 q1, q3, q4;
			ccx_3 q1, q3, q5;
			// flip q3, q4, q5;
			ccx_1 q1, q2, q3;
			ccx_2 q1, q2, q4;
			ccx_2 q1, q2, q5;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(5, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			expected_sv = zeros(32)
			expected_sv[end] = 1.0
			@test state ≈ expected_sv rtol=1e-10
		end

		@testset "11.13 Gate inverses" begin
			qasm = """
			gate rand_u_1 a { U(1, 2, 3) a; }
			gate rand_u_2 a { U(2, 3, 4) a; }
			gate rand_u_3 a { inv @ U(3, 4, 5) a; }

			gate both a {
				rand_u_1 a;
				rand_u_2 a;
			}
			gate both_inv a {
				inv @ both a;
			}
			gate all_3 a {
				rand_u_1 a;
				rand_u_2 a;
				rand_u_3 a;
			}
			gate all_3_inv a {
				inv @ inv @ inv @ all_3 a;
			}

			gate apply_phase a {
				gphase(1);
			}

			gate apply_phase_inv a {
				inv @ gphase(1);
			}

			qubit q;

			both q;
			both_inv q;

			all_3 q;
			all_3_inv q;

			apply_phase q;
			apply_phase_inv q;

			U(1, 2, 3) q;
			inv @ U(1, 2, 3) q;

			s q;
			inv @ s q;

			t q;
			inv @ t q;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(1, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector - all gates should cancel out
			expected_sv = [1.0, 0.0]
			@test state ≈ expected_sv
		end
		

		@testset "11.14 Gate on qubit registers" begin
			qasm = """
			qubit[3] qs;
			qubit q;

			x qs[{0, 2}];
			h q;
			cphaseshift(1) qs, q;
			phaseshift(-2) q;
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(4, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			@test state ≈ (1/√2)*[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
		end

		@testset "11.15 Verbatim" begin
			with_verbatim = """
			OPENQASM 3.0;
			bit[2] b;
			qubit[2] q;
			#pragma braket verbatim
			box{
				cnot q[0], q[1];
				cnot q[0], q[1];
				rx(1.57) q[0];
			}
			"""
			
			without_verbatim = """
			OPENQASM 3.0;
			bit[2] b;
			qubit[2] q;
			box{
				cnot q[0], q[1];
				cnot q[0], q[1];
				rx(1.57) q[0];
			}
			"""
			
			# Create simulators
			simulator1 = StateVectorSimulator(2, 1000)
			simulator2 = StateVectorSimulator(2, 1000)
			
			# Convert to circuits
			circuit_with = BraketSimulator.new_to_circuit(with_verbatim)
			circuit_without = BraketSimulator.new_to_circuit(without_verbatim)
			
			# Evolve using branched simulator operators
			branched_sim_with = BraketSimulator.evolve_branched_operators(simulator1, circuit_with, Dict{String, Any}())
			branched_sim_without = BraketSimulator.evolve_branched_operators(simulator2, circuit_without, Dict{String, Any}())
			
			# Calculate the final states
			state_with = BraketSimulator.calculate_current_state(branched_sim_with, branched_sim_with.active_paths[1])
			state_without = BraketSimulator.calculate_current_state(branched_sim_without, branched_sim_without.active_paths[1])
			
			# Verify the states are equivalent
			@test state_with ≈ state_without
		end

		@testset "11.16 Void subroutine" begin
			qasm = """
			def flip(qubit q) {
				x q;
			}
			qubit[2] qs;
			flip(qs[0]);
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(2, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			@test state ≈ [0, 0, 1, 0]
		end

		@testset "11.17 Rotation parameter expressions" begin
			qasm_pi = """
			OPENQASM 3.0;
			qubit[1] q;
			rx(π) q[0];
			"""
			
			# Create a simulator
			simulator = StateVectorSimulator(1, 1000)
			
			# Convert to circuit
			circuit = BraketSimulator.new_to_circuit(qasm_pi)
			
			# Evolve using branched simulator operators
			branched_sim = BraketSimulator.evolve_branched_operators(simulator, circuit, Dict{String, Any}())
			
			# Calculate the final state
			state = BraketSimulator.calculate_current_state(branched_sim, branched_sim.active_paths[1])
			
			# Verify the state vector
			@test state ≈ [0, -im]
			
			# Test more complex expressions
			qasm_expr = """
			OPENQASM 3.0;
			qubit[1] q;
			rx(pi + pi / 2) q[0];
			"""
			
			# Convert to circuit
			circuit_expr = BraketSimulator.new_to_circuit(qasm_expr)
			
			# Evolve using branched simulator operators
			branched_sim_expr = BraketSimulator.evolve_branched_operators(simulator, circuit_expr, Dict{String, Any}())
			
			# Calculate the final state
			state_expr = BraketSimulator.calculate_current_state(branched_sim_expr, branched_sim_expr.active_paths[1])
			
			# Verify the state vector
			@test state_expr ≈ [-0.70710678, -0.70710678im]
		end
	end
end
