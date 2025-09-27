using Test, LinearAlgebra, Random
using BraketSimulator

@testset "Branched Simulator Real Tests" begin
    @testset "Constructor and basic properties" begin
        base_sim = StateVectorSimulator(2, 1000)
        branched_sim = BranchedSimulator(base_sim)
        
        @test qubit_count(branched_sim) == 2
        @test device_id(branched_sim) == device_id(base_sim)
        @test branched_sim.shots == 1000
        @test length(branched_sim.active_paths) == 1
        @test length(branched_sim.states) == 1
        @test branched_sim.weights == [1.0]
        @test isempty(branched_sim.measurements[1])
    end
    
    @testset "Dynamic probability threshold" begin
        # Test different shot counts
        @test calculate_probability_threshold(0) == 0.001
        @test calculate_probability_threshold(50) == 0.05
        @test calculate_probability_threshold(500) == 0.01
        @test calculate_probability_threshold(5000) == 0.001
    end
    
    @testset "Path pruning" begin
        base_sim = StateVectorSimulator(2, 1000)
        branched_sim = BranchedSimulator(base_sim; max_paths=5)
        
        # Create more paths than max_paths
        for i in 1:10
            push!(branched_sim.states, copy(branched_sim.states[1]))
            push!(branched_sim.weights, 1.0 / (i + 1))
            push!(branched_sim.active_paths, length(branched_sim.states))
            push!(branched_sim.measurements, Dict{Int, Int}())
            push!(branched_sim.variables, Dict{String, Any}())
            push!(branched_sim.removed_qubits, Dict{Int, Int}())
            push!(branched_sim.qubit_mapping, Dict{Int, Int}())
        end
        
        # Should have 11 paths now (original + 10 new)
        @test length(branched_sim.active_paths) == 11
        
        # Prune paths
        BraketSimulator.prune_paths!(branched_sim)
        
        # Should have max_paths paths now
        @test length(branched_sim.active_paths) == branched_sim.max_paths
        
        # Weights should be normalized
        @test sum(branched_sim.weights[path_idx] for path_idx in branched_sim.active_paths) ≈ 1.0
    end
    
    @testset "Shot allocation" begin
        base_sim = StateVectorSimulator(2, 1000)
        branched_sim = BranchedSimulator(base_sim)
        
        # Create multiple paths with different weights
        push!(branched_sim.states, copy(branched_sim.states[1]))
        push!(branched_sim.weights, 0.3)
        push!(branched_sim.active_paths, 2)
        push!(branched_sim.measurements, Dict{Int, Int}())
        push!(branched_sim.variables, Dict{String, Any}())
        push!(branched_sim.removed_qubits, Dict{Int, Int}())
        push!(branched_sim.qubit_mapping, Dict{Int, Int}())
        
        push!(branched_sim.states, copy(branched_sim.states[1]))
        push!(branched_sim.weights, 0.2)
        push!(branched_sim.active_paths, 3)
        push!(branched_sim.measurements, Dict{Int, Int}())
        push!(branched_sim.variables, Dict{String, Any}())
        push!(branched_sim.removed_qubits, Dict{Int, Int}())
        push!(branched_sim.qubit_mapping, Dict{Int, Int}())
        
        # Adjust weight of first path
        branched_sim.weights[1] = 0.5
        
        # Allocate shots
        shot_allocation = BraketSimulator.allocate_shots(branched_sim)
        
        # Check that all shots are allocated
        @test sum(values(shot_allocation)) == 1000
        
        # Check that allocation is proportional to weights
        @test shot_allocation[1] ≈ 500 atol=5
        @test shot_allocation[2] ≈ 300 atol=5
        @test shot_allocation[3] ≈ 200 atol=5
    end
    
    @testset "Measurement and branching" begin
        # Create a simulator with a known state
        base_sim = StateVectorSimulator(1, 1000)
        # Apply Hadamard to create superposition
        BraketSimulator.apply_gate!(H(), base_sim.state_vector, 0)
        
        branched_sim = BranchedSimulator(base_sim)
        
        # Measure qubit 0
        BraketSimulator.measure_qubit!(branched_sim, 0)
        
        # Should have two paths now
        @test length(branched_sim.active_paths) == 2
        @test length(branched_sim.states) == 2
        
        # Check that measurements are recorded
        path1_idx = branched_sim.active_paths[1]
        path2_idx = branched_sim.active_paths[2]
        
        # One path should have measured 0, the other 1
        @test (branched_sim.measurements[path1_idx][0] == 0 && branched_sim.measurements[path2_idx][0] == 1) ||
              (branched_sim.measurements[path1_idx][0] == 1 && branched_sim.measurements[path2_idx][0] == 0)
        
        # Weights should be proportional to probabilities
        @test branched_sim.weights[path1_idx] ≈ 0.5
        @test branched_sim.weights[path2_idx] ≈ 0.5
    end
    
    @testset "Classical variable handling" begin
        base_sim = StateVectorSimulator(1, 1000)
        branched_sim = BranchedSimulator(base_sim)
        
        # Set a variable
        BraketSimulator.set_variable!(branched_sim, 1, "test_var", 42)
        
        # Get the variable
        @test BraketSimulator.get_variable(branched_sim, 1, "test_var") == 42
        
        # Get a non-existent variable
        @test BraketSimulator.get_variable(branched_sim, 1, "non_existent") === nothing
        @test BraketSimulator.get_variable(branched_sim, 1, "non_existent", "default") == "default"
    end
    
    @testset "Memory optimization with statevector reduction" begin
        # Create a simulator with 3 qubits
        base_sim = StateVectorSimulator(3, 1000)
        branched_sim = BranchedSimulator(base_sim)
        
        # Initial state should have 2^3 = 8 amplitudes
        @test length(branched_sim.states[1]) == 8
        
        # Measure qubit 0
        BraketSimulator.measure_qubit!(branched_sim, 0)
        
        # Should have two paths now
        @test length(branched_sim.active_paths) == 2
        
        # Each path's statevector should be reduced to 2^2 = 4 amplitudes
        for path_idx in branched_sim.active_paths
            @test length(branched_sim.states[path_idx]) == 4
            
            # Check that qubit 0 is recorded as removed
            @test haskey(branched_sim.removed_qubits[path_idx], 0)
            
            # Check that the qubit mapping is updated
            @test !haskey(branched_sim.qubit_mapping[path_idx], 0)
            @test haskey(branched_sim.qubit_mapping[path_idx], 1)
            @test haskey(branched_sim.qubit_mapping[path_idx], 2)
        end
        
        # Measure qubit 1
        BraketSimulator.measure_qubit!(branched_sim, 1)
        
        # Should have four paths now
        @test length(branched_sim.active_paths) == 4
        
        # Each path's statevector should be reduced to 2^1 = 2 amplitudes
        for path_idx in branched_sim.active_paths
            @test length(branched_sim.states[path_idx]) == 2
            
            # Check that qubits 0 and 1 are recorded as removed
            @test haskey(branched_sim.removed_qubits[path_idx], 0)
            @test haskey(branched_sim.removed_qubits[path_idx], 1)
            
            # Check that the qubit mapping is updated
            @test !haskey(branched_sim.qubit_mapping[path_idx], 0)
            @test !haskey(branched_sim.qubit_mapping[path_idx], 1)
            @test haskey(branched_sim.qubit_mapping[path_idx], 2)
        end
    end
    
    @testset "Statevector expansion when operating on measured qubits" begin
        # Create a simulator with 2 qubits
        base_sim = StateVectorSimulator(2, 1000)
        branched_sim = BranchedSimulator(base_sim)
        
        # Measure qubit 0
        BraketSimulator.measure_qubit!(branched_sim, 0)
        
        # Should have two paths now
        @test length(branched_sim.active_paths) == 2
        
        # Each path's statevector should be reduced to 2^1 = 2 amplitudes
        for path_idx in branched_sim.active_paths
            @test length(branched_sim.states[path_idx]) == 2
            @test haskey(branched_sim.removed_qubits[path_idx], 0)
        end
        
        # Now expand the statevector to reincorporate qubit 0
        for path_idx in branched_sim.active_paths
            BraketSimulator.expand_statevector!(branched_sim, path_idx, 0)
        end
        
        # Each path's statevector should be expanded back to 2^2 = 4 amplitudes
        for path_idx in branched_sim.active_paths
            @test length(branched_sim.states[path_idx]) == 4
            @test !haskey(branched_sim.removed_qubits[path_idx], 0)
            @test haskey(branched_sim.qubit_mapping[path_idx], 0)
        end
    end
    
    @testset "Mid-circuit measurement with OpenQASM" begin
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
        
        # Parse the program
        program = BraketSimulator.new_to_circuit(qasm_source)
        
        # Create a simulator
        simulator = StateVectorSimulator(2, 1000)
        
        # Create a branched simulator and evolve the program
        branched_sim = BraketSimulator.evolve_branched!(simulator, program, Dict{String, Any}())
        
        # Should have two paths now
        @test length(branched_sim.active_paths) == 2
        
        # Check that measurements are recorded
        path1_idx = branched_sim.active_paths[1]
        path2_idx = branched_sim.active_paths[2]
        
        # One path should have measured 0, the other 1
        @test (branched_sim.measurements[path1_idx][0] == 0 && branched_sim.measurements[path2_idx][0] == 1) ||
              (branched_sim.measurements[path1_idx][0] == 1 && branched_sim.measurements[path2_idx][0] == 0)
        
        # Allocate shots
        shot_allocation = BraketSimulator.allocate_shots(branched_sim)
        
        # Check that all shots are allocated
        @test sum(values(shot_allocation)) == 1000
        
        # Check that allocation is roughly equal (since H gives 50/50)
        @test shot_allocation[path1_idx] ≈ 500 atol=50
        @test shot_allocation[path2_idx] ≈ 500 atol=50
    end
end
