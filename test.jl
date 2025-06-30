using BraketSimulator

# Test dynamic qubit allocation in the branched simulator
function test_dynamic_qubit_allocation()
	println("Testing dynamic qubit allocation...")

	# Create a simulator with 100 shots
	# The 2 initial qubits is irrelevant as the branched simulator assigns qubits dynamically
	simulator = StateVectorSimulator(2, 100)

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

	# Evolve the program using the branched simulator
	branched_sim = BraketSimulator.evolve_branched_operators(simulator, BraketSimulator.new_to_circuit(qasm_source), Dict{String, Any}())
	println(branched_sim.instruction_sequences)
	println(branched_sim.measurements)
	println(branched_sim.variables)
	println(branched_sim.active_paths)
	paths_by_count = Dict{Int, Vector{Int}}()

	for path_idx in branched_sim.active_paths
		count = BraketSimulator.get_variable(branched_sim, path_idx, "count").val
		if !haskey(paths_by_count, count)
			paths_by_count[count] = Int[]
		end
		push!(paths_by_count[count], path_idx)
	end


	println(paths_by_count)
    println(BraketSimulator.calculate_current_state(branched_sim, 3))
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
