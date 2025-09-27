using Test, LinearAlgebra
using BraketSimulator: apply_projection, expand_state, get_measurement_probabilities, remove_bit

@testset "Standalone Measurement Operations" begin
	@testset "StateVector Operations" begin
		@testset "Measurement Probabilities" begin
			# Test |0⟩ state
			sv = zeros(ComplexF64, 2)
			sv[1] = 1.0
			probs = get_measurement_probabilities(sv, 0)
			@test probs ≈ [1.0, 0.0]

			# Test |1⟩ state
			sv = zeros(ComplexF64, 2)
			sv[2] = 1.0
			probs = get_measurement_probabilities(sv, 0)
			@test probs ≈ [0.0, 1.0]

			# Test |+⟩ state
			sv = ones(ComplexF64, 2) ./ sqrt(2)
			probs = get_measurement_probabilities(sv, 0)
			@test probs ≈ [0.5, 0.5]

			# Test 2-qubit |00⟩ state
			sv = zeros(ComplexF64, 4)
			sv[1] = 1.0
			probs_q0 = get_measurement_probabilities(sv, 0)
			probs_q1 = get_measurement_probabilities(sv, 1)
			@test probs_q0 ≈ [1.0, 0.0]
			@test probs_q1 ≈ [1.0, 0.0]

			# Test 2-qubit Bell state (|00⟩ + |11⟩)/√2
			sv = zeros(ComplexF64, 4)
			sv[1] = 1.0/sqrt(2)
			sv[4] = 1.0/sqrt(2)
			probs_q0 = get_measurement_probabilities(sv, 0)
			probs_q1 = get_measurement_probabilities(sv, 1)
			@test probs_q0 ≈ [0.5, 0.5]
			@test probs_q1 ≈ [0.5, 0.5]
		end

		@testset "Projection" begin
			# Test 2-qubit Bell state (|00⟩ + |11⟩)/√2 projection
			sv = zeros(ComplexF64, 4)
			sv[1] = 1.0/sqrt(2)
			sv[4] = 1.0/sqrt(2)

			# Project first qubit to |0⟩
			sv_copy = copy(sv)
			sv_copy = apply_projection(sv_copy, 0, 0)
			@test sv_copy ≈ [1.0, 0.0]

			# Project first qubit to |1⟩
			sv_copy = copy(sv)
			sv_copy = apply_projection(sv_copy, 0, 1)
			@test sv_copy ≈ [0.0, 1.0]
		end

		@testset "Expand State" begin
			# Test expanding 1-qubit |0⟩ state by adding qubit 0 in |0⟩ state
			sv = zeros(ComplexF64, 2)
			sv[1] = 1.0
			new_sv = expand_state(sv, 0, 0)
			@test length(new_sv) == 4
			@test new_sv ≈ [1.0, 0.0, 0.0, 0.0]

			# Test expanding 1-qubit |1⟩ state by adding qubit 0 in |1⟩ state
			sv = zeros(ComplexF64, 2)
			sv[2] = 1.0
			new_sv = expand_state(sv, 0, 1)
			@test length(new_sv) == 4
			@test new_sv ≈ [0.0, 0.0, 0.0, 1.0]

            # Test expanding 2-qubit |00⟩ state by adding qubit 0 in |1⟩ state
			sv = zeros(ComplexF64, 4)
			sv[1] = 1.0
			new_sv = expand_state(sv, 0, 1)
			@test length(new_sv) == 8
			@test new_sv ≈ [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]

            # Test expanding 2-qubit |00⟩ state by adding qubit 1 in |1⟩ state
			sv = zeros(ComplexF64, 4)
			sv[1] = 1.0
			new_sv = expand_state(sv, 1, 1)
			@test length(new_sv) == 8
			@test new_sv ≈ [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
		end
	end

	@testset "DensityMatrix Operations" begin
		@testset "Measurement Probabilities" begin
			# Test |0⟩⟨0| state
			dm = zeros(ComplexF64, 2, 2)
			dm[1, 1] = 1.0
			probs = get_measurement_probabilities(dm, 0)
			@test probs ≈ [1.0, 0.0]

			# Test |1⟩⟨1| state
			dm = zeros(ComplexF64, 2, 2)
			dm[2, 2] = 1.0
			probs = get_measurement_probabilities(dm, 0)
			@test probs ≈ [0.0, 1.0]

			# Test maximally mixed state
			dm = Matrix{ComplexF64}(I, 2, 2) ./ 2
			probs = get_measurement_probabilities(dm, 0)
			@test probs ≈ [0.5, 0.5]
		end

		@testset "Projection" begin
			# Test 2-qubit Bell state density matrix projection
			bell_sv = zeros(ComplexF64, 4)
			bell_sv[1] = 1.0/sqrt(2)
			bell_sv[4] = 1.0/sqrt(2)
			bell_dm = bell_sv * bell_sv'

			# Project first qubit to |0⟩
			dm_copy = copy(bell_dm)
			dm_copy = apply_projection(dm_copy, 0, 0)
			# Should be a 2x2 matrix representing |0⟩⟨0| for the remaining qubit
			@test size(dm_copy) == (2, 2)
			@test dm_copy ≈ [1.0 0.0; 0.0 0.0]

			# Project first qubit to |1⟩
			dm_copy = copy(bell_dm)
			dm_copy = apply_projection(dm_copy, 0, 1)
			# Should be a 2x2 matrix representing |1⟩⟨1| for the remaining qubit
			@test size(dm_copy) == (2, 2)
			@test dm_copy ≈ [0.0 0.0; 0.0 1.0]

			# Test 3-qubit GHZ state projection
			ghz_sv = zeros(ComplexF64, 8)
			ghz_sv[1] = 1.0/sqrt(2)  # |000⟩
			ghz_sv[8] = 1.0/sqrt(2)  # |111⟩
			ghz_dm = ghz_sv * ghz_sv'

			# Project first qubit to |0⟩
			dm_copy = copy(ghz_dm)
			dm_copy = apply_projection(dm_copy, 0, 0)
			# Should be a 4x4 matrix representing |00⟩⟨00| for the remaining qubits
			@test size(dm_copy) == (4, 4)
			expected = zeros(ComplexF64, 4, 4)
			expected[1, 1] = 1.0
			@test dm_copy ≈ expected

			# Project middle qubit to |1⟩
			dm_copy = copy(ghz_dm)
			dm_copy = apply_projection(dm_copy, 1, 1)
			# Should be a 4x4 matrix representing a mixture of |01⟩⟨01| and |10⟩⟨10|
			@test size(dm_copy) == (4, 4)
			expected = zeros(ComplexF64, 4, 4)
			expected[4, 4] = 1.0
			@test dm_copy ≈ expected
		end

		@testset "Integration Tests" begin
			# Test workflow: create Bell state, project (which now includes reduction), expand
			bell_sv = zeros(ComplexF64, 4)
			bell_sv[1] = 1.0/sqrt(2)
			bell_sv[4] = 1.0/sqrt(2)
			bell_dm = bell_sv * bell_sv'

			# Project first qubit to |0⟩ (now includes reduction)
			dm_copy = copy(bell_dm)
			dm_copy = apply_projection(dm_copy, 0, 0)
			@test size(dm_copy) == (2, 2)
			@test dm_copy ≈ [1.0 0.0; 0.0 0.0]

			# Expand state to reincorporate qubit 0
			expanded_dm = expand_state(copy(dm_copy), 0, 0)
			@test size(expanded_dm) == (4, 4)
			expected = zeros(ComplexF64, 4, 4)
			expected[1, 1] = 1.0
			@test expanded_dm ≈ expected

			# Test with a 3-qubit GHZ state
			ghz_sv = zeros(ComplexF64, 8)
			ghz_sv[1] = 1.0/sqrt(2)  # |000⟩
			ghz_sv[8] = 1.0/sqrt(2)  # |111⟩
			ghz_dm = ghz_sv * ghz_sv'

			# Project and reduce qubit 0 to |0⟩
			dm_copy = copy(ghz_dm)
			dm_copy = apply_projection(dm_copy, 0, 0)
			@test size(dm_copy) == (4, 4)
			expected_reduced = zeros(ComplexF64, 4, 4)
			expected_reduced[1, 1] = 1.0
			@test dm_copy ≈ expected_reduced

			# Project and reduce qubit 1 to |1⟩ from the already reduced state
			further_reduced = apply_projection(dm_copy, 1, 0)
			@test size(further_reduced) == (2, 2)
			expected_further_reduced = zeros(ComplexF64, 2, 2)
			expected_further_reduced[1, 1] = 1.0
			@test further_reduced ≈ expected_further_reduced

			# Test sequential measurements with expansion
			# Start with GHZ state, measure qubit 0 to |0⟩, then expand to add qubit 0 back as |1⟩
			dm_copy = copy(ghz_dm)
			dm_copy = apply_projection(dm_copy, 0, 0)
			expanded_dm = expand_state(copy(dm_copy), 0, 1)
			@test size(expanded_dm) == (8, 8)
			
			# The expanded state should be |100⟩⟨100| (since we measured qubit 0 to |0⟩ but reinserted it as |1⟩)
			expected_expanded = zeros(ComplexF64, 8, 8)
			expected_expanded[5, 5] = 1.0
			@test expanded_dm ≈ expected_expanded
		end
	end
end
