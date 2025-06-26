"""
    DensityMatrixSimulator{T, S<:AbstractMatrix{T}} <: AbstractSimulator

Simulator representing evolution of a density matrix of type `S`, with element type `T`. Density matrix simulators should be used to simulate circuits with noise.
"""
mutable struct DensityMatrixSimulator{T,S} <:
               AbstractSimulator where {T,S<:AbstractDensityMatrix{T}}
    density_matrix::S
    qubit_count::Int
    shots::Int
    shot_buffer::Vector{Int}
    _alias::Vector{Int}
    _ap::Vector{Float64}
    _larges::Vector{Int}
    _smalls::Vector{Int}
    _density_matrix_after_observables::S
    function DensityMatrixSimulator{T,S}(
        density_matrix::S,
        qubit_count::Int,
        shots::Int,
    ) where {T,S<:AbstractDensityMatrix{T}}
        shot_buffer = Vector{Int}(undef, shots)
        ap_len  = ap_size(shots, qubit_count)
        _ap     = zeros(Float64, ap_len)
        _alias  = zeros(Int, ap_len)
        _larges = zeros(Int, ap_len)
        _smalls = zeros(Int, ap_len)
        return new(
            density_matrix,
            qubit_count,
            shots,
            shot_buffer,
            _alias,
            _ap,
            _larges,
            _smalls,
            S(undef, 0, 0),
        )
    end
end
function init(
    t::Type{S},
    qubit_count::Int,
) where {T<:Complex,S<:AbstractMatrix{T}}
    dm = t(undef, 2^qubit_count, 2^qubit_count)
    fill!(dm, zero(T))
    dm[1, 1] = one(T)
    return dm
end
function DensityMatrixSimulator{T,S}(
    qubit_count::Int,
    shots::Int,
) where {T,S<:AbstractDensityMatrix{T}}
    dm = init(S, qubit_count)
    return DensityMatrixSimulator{T,S}(dm, qubit_count, shots)
end
"""
    DensityMatrixSimulator([::T], qubit_count::Int, shots::Int) -> DensityMatrixSimulator{T, Matrix{T}}

Create a `DensityMatrixSimulator` with `2^qubit_count x 2^qubit_count` elements and `shots` shots to be measured. The default element type is `ComplexF64`.
"""
DensityMatrixSimulator(::Type{T}, qubit_count::Int, shots::Int) where {T<:Number} =
    DensityMatrixSimulator{T,Matrix{T}}(qubit_count, shots)
DensityMatrixSimulator(qubit_count::Int, shots::Int) =
    DensityMatrixSimulator(ComplexF64, qubit_count, shots)
qubit_count(dms::DensityMatrixSimulator) = dms.qubit_count
"""
    properties(svs::DensityMatrixSimulator) -> GateModelSimulatorDeviceCapabilities

Query the properties and capabilities of a `DensityMatrixSimulator`, including which gates and result types are supported and the minimum and maximum shot and qubit counts.
"""
properties(d::DensityMatrixSimulator) = dm_props
supported_operations(d::DensityMatrixSimulator, ::Val{:OpenQASM}) = dm_props.action["braket.ir.openqasm.program"].supportedOperations
supported_operations(d::DensityMatrixSimulator, ::Val{:JAQCD}) = dm_props.action["braket.ir.jaqcd.program"].supportedOperations
supported_operations(d::DensityMatrixSimulator) = supported_operations(d::DensityMatrixSimulator, Val(:OpenQASM))
supported_result_types(d::DensityMatrixSimulator, ::Val{:OpenQASM}) = dm_props.action["braket.ir.openqasm.program"].supportedResultTypes
supported_result_types(d::DensityMatrixSimulator, ::Val{:JAQCD}) = dm_props.action["braket.ir.jaqcd.program"].supportedResultTypes
supported_result_types(d::DensityMatrixSimulator) = supported_result_types(d::DensityMatrixSimulator, Val(:OpenQASM))
device_id(dms::DensityMatrixSimulator) = "braket_dm_v2"
name(dms::DensityMatrixSimulator) = "DensityMatrixSimulator"
Base.show(io::IO, dms::DensityMatrixSimulator) =
    print(io, "DensityMatrixSimulator(qubit_count=$(qubit_count(dms)), shots=$(dms.shots))")
Base.similar(
    dms::DensityMatrixSimulator{T,S};
    shots::Int = dms.shots,
) where {T,S<:AbstractDensityMatrix{T}} =
    DensityMatrixSimulator{T,S}(dms.qubit_count, shots)
Base.copy(dms::DensityMatrixSimulator{T,S}) where {T,S<:AbstractDensityMatrix{T}} =
    DensityMatrixSimulator{T,S}(deepcopy(dms.density_matrix), dms.qubit_count, dms.shots)
function Base.copyto!(
    dst::DensityMatrixSimulator{T,S},
    src::DensityMatrixSimulator{T,S},
) where {T,S}
    copyto!(dst.density_matrix, src.density_matrix)
    return dst
end

function reinit!(
    dms::DensityMatrixSimulator{T,S},
    qubit_count::Int,
    shots::Int,
) where {T,S<:AbstractDensityMatrix{T}}
    n_amps = 2^qubit_count
    if size(dms.density_matrix) != (n_amps, n_amps)
        dms.density_matrix = S(undef, n_amps, n_amps)
        ap_len = ap_size(shots, qubit_count)
        resize!(dms._alias, ap_len)
        resize!(dms._ap, ap_len)
        resize!(dms._larges, ap_len)
        resize!(dms._smalls, ap_len)
    end
    if dms.shots != shots
        ap_len = ap_size(shots, qubit_count)
        resize!(dms._alias, ap_len)
        resize!(dms._ap, ap_len)
        resize!(dms._larges, ap_len)
        resize!(dms._smalls, ap_len)
        resize!(dms.shot_buffer, shots)
    end
    fill!(dms.density_matrix, zero(T))
    dms._ap        .= zero(Float64)
    dms._alias     .= zero(Int)
    dms._larges    .= zero(Int)
    dms._smalls    .= zero(Int)
    dms.qubit_count = qubit_count
    dms.shots       = shots
    dms.density_matrix[1, 1] = one(T)
    dms._density_matrix_after_observables = S(undef, 0, 0)
    return
end

function _evolve_op!(
    dms::DensityMatrixSimulator{T,S},
    op::G,
    target::Int...,
) where {T<:Complex,S<:AbstractDensityMatrix{T},G<:Gate}
    reshaped_dm = reshape(dms.density_matrix, length(dms.density_matrix))
    apply_gate!(Val(false), op, reshaped_dm, target...)
    # applies the *conjugate* of `op` to the "qubits" corresponding to the *transpose*
    apply_gate!(Val(true),  op, reshaped_dm, (dms.qubit_count .+ target)...)
    return
end

function _evolve_op!(
    dms::DensityMatrixSimulator{T,S},
    op::N,
    target::Int...,
) where {T<:Complex,S<:AbstractDensityMatrix{T},N<:Noise}
    apply_noise!(op, dms.density_matrix, target...)
end
# Measure, barrier, reset, and delay operators are no-ops for now as
# measurement is handled at the end of simulation, in the results
# computation step. If/when mid-circuit # measurement is supported,
# this operation will collapse the density # matrix on the measured qubits.
# Barrier, reset, and delay are also to-do implementations.
_evolve_op!(dms::DensityMatrixSimulator{T,S}, m::Measure, args...) where {T<:Complex,S<:AbstractDensityMatrix{T}} = return
_evolve_op!(dms::DensityMatrixSimulator{T,S}, b::Barrier, args...) where {T<:Complex,S<:AbstractDensityMatrix{T}} = return
_evolve_op!(dms::DensityMatrixSimulator{T,S}, r::Reset, args...) where {T<:Complex,S<:AbstractDensityMatrix{T}} = return
_evolve_op!(dms::DensityMatrixSimulator{T,S}, d::Delay, args...) where {T<:Complex,S<:AbstractDensityMatrix{T}} = return

"""
    evolve!(dms::DensityMatrixSimulator{T, S<:AbstractMatrix{T}}, operations::Vector{Instruction}) -> DensityMatrixSimulator{T, S}

Apply each operation of `operations` in-place to the density matrix contained in `dms`.

Effectively, perform the operation:

`` \\hat{\\rho} \\to \\hat{A}^\\dag \\hat{\\rho} \\hat{A} ``

for each operation ``\\hat{A}`` in `operations`.
"""
function evolve!(
    dms::DensityMatrixSimulator{T,S},
    operations,
)::DensityMatrixSimulator{T,S} where {T<:Complex,S<:AbstractDensityMatrix{T}}
    for operation in operations
        # use this to dispatch on Gates vs Noises
        _evolve_op!(dms, operation.operator, operation.target...)
    end
    return dms
end

function apply_observable!(
    gate::G,
    dm::S,
    targets,
) where {T<:Complex,S<:AbstractDensityMatrix{T}, G<:Gate}
    reshaped_dm = reshape(dm, length(dm))
    foreach(target->apply_gate!(gate, reshaped_dm, target), targets)
    return dm
end
function apply_observable!(
    observable::Observables.HermitianObservable,
    dm::AbstractDensityMatrix{T},
    targets::Int...,
) where {T<:Complex}
    n_amps    = size(dm, 1)
    n_qubits  = Int(log2(n_amps))
    endian_ts = n_qubits - 1 .- collect(targets)
    n_targets = length(endian_ts)
    o_mat     = transpose(observable.matrix)

    ordered_ts = sort(endian_ts)
    flip_list = map(0:2^n_targets-1) do target_amp
        to_flip = Bool[(((1 << target_ix) & target_amp) >> target_ix) for target_ix = 0:n_targets-1]
        return ordered_ts[to_flip]
    end
    slim_size = div(n_amps, 2^n_targets)
    Threads.@threads for raw_ix = 0:(slim_size^2)-1
        ix = div(raw_ix, slim_size)
        jx = mod(raw_ix, slim_size)
        padded_ix = pad_bits(ix, ordered_ts)
        padded_jx = pad_bits(jx, ordered_ts)
        ixs = map(flip_list) do bits_to_flip
            flip_bits(padded_ix, bits_to_flip) + 1
        end
        jxs = map(flip_list) do bits_to_flip 
            flip_bits(padded_jx, bits_to_flip) + 1
        end
        @views @inbounds begin
            elems = dm[jxs[:], ixs[:]]
            dm[jxs[:], ixs[:]] = o_mat * elems
        end
    end
    return dm
end

function state_with_observables(dms::DensityMatrixSimulator)
    isempty(dms._density_matrix_after_observables) &&
        error("observables have not been applied.")
    return dms._density_matrix_after_observables
end

function apply_observables!(dms::DensityMatrixSimulator, observables)
    !isempty(dms._density_matrix_after_observables) &&
        error("observables have already been applied.")
    operations = mapreduce(obs->diagonalizing_gates(obs...), vcat, observables)
    dms._density_matrix_after_observables = deepcopy(dms.density_matrix)
    reshaped_dm = reshape(dms._density_matrix_after_observables, length(dms.density_matrix))
    for operation in operations
        apply_gate!(Val(false), operation.operator, reshaped_dm, operation.target...)
        apply_gate!(Val(true), operation.operator, reshaped_dm, (dms.qubit_count .+ operation.target)...)
    end
    return dms
end

"""
    expectation(dms::DensityMatrixSimulator, observable::Observables.Observable, targets::Int...) -> Float64

Compute the exact (`shots=0`) expectation value of `observable` applied to `targets`
given the evolved density matrix in `dms`. In other words, compute

``\\mathrm{Tr}\\left(\\hat{O}\\hat{\\rho}\\right)``.
"""
function expectation(
    dms::DensityMatrixSimulator,
    observable::Observables.Observable,
    targets::Int...,
)
    dm_copy = apply_observable(observable, dms.density_matrix, targets...)
    return real(sum(diag(dm_copy)))
end
density_matrix(dms::DensityMatrixSimulator) = copy(dms.density_matrix)
probabilities(dms::DensityMatrixSimulator)  = real.(diag(dms.density_matrix))

function swap_bits(ix::Int, qubit_map::Dict{Int,Int})
    # only flip 01 and 10
    # flipping 00 and 11 will not change the final
    # integer index
    for (in_q, out_q) in qubit_map
        if in_q < out_q # avoid flipping twice
            in_val  = ((1 << in_q) & ix) >> in_q
            out_val = ((1 << out_q) & ix) >> out_q
            if in_val != out_val
                ix = flip_bits(ix, (in_q, out_q))
            end
        end
    end
    return ix
end

function partial_trace(
    ρ::AbstractMatrix{ComplexF64},
    output_qubits = collect(0:Int(log2(size(ρ, 1)))-1),
)
    isempty(output_qubits) && return sum(diag(ρ))
    n_amps   = size(ρ, 1)
    n_qubits = Int(log2(n_amps))
    length(unique(output_qubits)) == n_qubits && return ρ

    qubits        = setdiff(collect(0:n_qubits-1), output_qubits)
    endian_qubits = sort(n_qubits .- qubits .- 1)
    qubit_combos  = vcat([Int[]], collect(combinations(endian_qubits)))
    final_ρ_dim   = 2^(n_qubits - length(qubits))
    final_ρ       = zeros(ComplexF64, final_ρ_dim, final_ρ_dim)
    # handle possibly permuted targets
    needs_perm     = !issorted(output_qubits)
    final_n_qubits = length(output_qubits)
    output_qubit_mapping = if needs_perm
            original_outputs = final_n_qubits .- output_qubits .- 1
            permuted_outputs = final_n_qubits .- collect(0:final_n_qubits-1) .- 1
            Dict(zip(original_outputs, permuted_outputs))
        else
            Dict{Int,Int}()
        end
    for raw_ix = 0:length(final_ρ)-1
        ix = div(raw_ix, size(final_ρ, 1))
        jx = mod(raw_ix, size(final_ρ, 1))
        padded_ix = pad_bits(ix, endian_qubits)
        padded_jx = pad_bits(jx, endian_qubits)
        flipped_inds = Vector{CartesianIndex{2}}(undef, length(qubit_combos))
        for (c_ix, flipped_qubits) in enumerate(qubit_combos)
            flipped_ix = flip_bits(padded_ix, flipped_qubits)
            flipped_jx = flip_bits(padded_jx, flipped_qubits)
            flipped_inds[c_ix] = CartesianIndex{2}(flipped_ix + 1, flipped_jx + 1)
        end
        # if the output qubits weren't in sorted order, we need to permute the
        # final indices of ρ to match the desired qubit mapping
        out_ix = needs_perm ? swap_bits(ix, output_qubit_mapping) : ix
        out_jx = needs_perm ? swap_bits(jx, output_qubit_mapping) : jx
        @views @inbounds final_ρ[out_ix+1, out_jx+1] = sum(ρ[flipped_inds])
    end
    return final_ρ
end
partial_trace(sim::DensityMatrixSimulator, output_qubits = collect(0:qubit_count(sim)-1)) = partial_trace(density_matrix(sim), output_qubits)

##################################
# Measurement Related Operations #
##################################


"""
    get_measurement_probabilities(state::AbstractMatrix{<:Complex}, qubit::Int)

Calculate measurement probabilities for a qubit in a density matrix.
Returns [P(0), P(1)] for the given qubit.
"""
function get_measurement_probabilities(state::AbstractMatrix{<:Complex}, qubit::Int)
    n_qubits = Int(log2(size(state, 1)))
    endian_qubit = n_qubits - qubit - 1
    
    # Calculate probabilities for 0 and 1 outcomes
    prob_0 = 0.0
    prob_1 = 0.0
    
    # For density matrix, we need to use the diagonal elements
    for i in 0:(size(state, 1)-1)
        # Check if the qubit is 0 or 1 in this basis state
        qubit_val = (i >> endian_qubit) & 1
        prob = real(state[i+1, i+1])  # Probability from diagonal element
        
        if qubit_val == 0
            prob_0 += prob
        else
            prob_1 += prob
        end
    end
    
    return [prob_0, prob_1]
end

"""
    get_measurement_probabilities(simulator::DensityMatrixSimulator, qubit::Int)

Calculate measurement probabilities for a qubit in a density matrix simulator.
Returns [P(0), P(1)] for the given qubit.
"""
function get_measurement_probabilities(simulator::DensityMatrixSimulator, qubit::Int)
    return get_measurement_probabilities(simulator.density_matrix, qubit)
end

"""
    apply_projection(state::AbstractMatrix{<:Complex}, qubit::Int, outcome::Int)

Apply measurement projection to a density matrix.
This collapses the state to the subspace corresponding to the measurement outcome
and returns the reduced density matrix. If all qubits have been measured, the returned matrix is empty.
Throws an error if the qubit outcome isn't possible
"""
function apply_projection(state::AbstractMatrix{<:Complex}, qubit::Int, outcome::Int)
    n_qubits = Int(log2(size(state, 1)))
    endian_qubit = n_qubits - qubit - 1
    dim = size(state, 1)

    # If the state is empty, then there is nothing to measure!
    if dim == 0
        throw("There are no qubits to be measured in the density matrix")
    end

    # If there is only one qubit present in the state, then it must be measured and the projection
    # removes the qubit from the state, so you are just left with an empty state vector
    if n_qubits == 1
        return Matrix{ComplexF64}(undef, 0, 0)
    end
    
    # Calculate normalization factor (trace of projected state)
    trace_val = get_measurement_probabilities(state, qubit)[outcome+1]
    
    # This means that the qubit outcome is impossible
    if trace_val <= 0
        throw("Qubit measurement outcome not possible, please check measurement probability before applying the projection!")
    end
    
    # Create reduced density matrix (half the size)
    reduced_dim = dim ÷ 2
    reduced_state = zeros(ComplexF64, reduced_dim, reduced_dim)
    
    # Fill the reduced density matrix
    for i in 0:(dim-1)
        i_qubit_val = (i >> endian_qubit) & 1
        if i_qubit_val != outcome
            continue  # Skip elements that don't match the outcome
        end
        
        # Remove the measured qubit from index i
        reduced_i = remove_bit(i, endian_qubit)
        
        for j in 0:(dim-1)
            j_qubit_val = (j >> endian_qubit) & 1
            if j_qubit_val != outcome
                continue  # Skip elements that don't match the outcome
            end
            
            # Remove the measured qubit from index j
            reduced_j = remove_bit(j, endian_qubit)
            
            # Copy and normalize the element
            reduced_state[reduced_i+1, reduced_j+1] = state[i+1, j+1] / trace_val
        end
    end
    
    return reduced_state
end

"""
    remove_bit(index::Int, position::Int) -> Int

Remove a bit at the specified position in the binary representation of index.
Used for density matrix index manipulation during reduction.
"""
function remove_bit(index::Int, position::Int)
    # Extract bits before the position
    mask_before = (1 << position) - 1
    bits_before = index & mask_before
    
    # Extract bits after the position, shifted right
    bits_after = (index >> (position + 1)) << position
    
    # Combine parts
    return bits_before | bits_after
end

"""
    _apply_reset(state::AbstractMatrix{<:Complex}, target::Int)

Reset a qubit to |0⟩ state in a density matrix.
This projects the state to the |0⟩ state of the target qubit.
"""
function _apply_reset(state::AbstractMatrix{<:Complex}, target::Int)
    # First measure the qubit (project to either |0⟩ or |1⟩)
    probs = get_measurement_probabilities(state, target)
    outcome = rand() < probs[1] ? 0 : 1
    apply_projection!(state, target, outcome)
    
    # If outcome was |1⟩, apply X gate to make it |0⟩
    if outcome == 1
        n_qubits = Int(log2(size(state, 1)))
        endian_target = n_qubits - target - 1
        dim = size(state, 1)
        
        # Create X gate matrix for the target qubit
        x_matrix = zeros(ComplexF64, dim, dim)
        for i in 0:(dim-1)
            bit_val = (i >> endian_target) & 1
            flipped_i = i ⊻ (1 << endian_target)
            x_matrix[flipped_i+1, i+1] = 1.0
        end
        
        # Apply X gate: ρ -> X ρ X†
        state .= x_matrix * state * adjoint(x_matrix)
    end
    
    return state
end

"""
    _apply_reset!(simulator::DensityMatrixSimulator, target::Int)

Reset a qubit to |0⟩ state in a density matrix simulator.
This projects the state to the |0⟩ state of the target qubit.
"""
function _apply_reset(simulator::DensityMatrixSimulator, target::Int)
    return _apply_reset(simulator.density_matrix, target)
end


"""
    expand_state(state::AbstractMatrix{<:Complex}, qubit::Int, outcome::Int)

Expand the density matrix to reincorporate a previously measured qubit.
This is needed when an operation is applied to a qubit that was previously measured, qubit reuse.
Returns the expanded density matrix
"""
function expand_state(state::AbstractMatrix{<:Complex}, qubit::Int, outcome::Int)
    n_qubits = Int(log2(size(state, 1)))
    
    # Create a new density matrix with expanded size
    new_dim = 2^(n_qubits + 1)
    new_state = zeros(ComplexF64, new_dim, new_dim)
    
    # Determine where to insert the qubit in the new density matrix
    endian_qubit = n_qubits - qubit
    
    # Map from original indices to expanded indices with the measured qubit inserted
    for i in 0:(size(state, 1)-1)
        for j in 0:(size(state, 1)-1)
            # Calculate new indices with the qubit inserted
            new_i = insert_bit(i, endian_qubit, outcome)
            new_j = insert_bit(j, endian_qubit, outcome)
            
            # Copy the matrix element
            new_state[new_i+1, new_j+1] = state[i+1, j+1]
        end
    end
    
    return new_state
end

"""
    insert_bit(index::Int, position::Int, bit::Int) -> Int

Insert a bit at the specified position in the binary representation of index.
Used for density matrix index manipulation during expansion/reduction.
"""
function insert_bit(index::Int, position::Int, bit::Int)
    # Create mask for bits before the position
    mask_before = (1 << position) - 1
    
    # Create mask for bits after the position
    mask_after = ~mask_before
    
    # Extract bits before and after the position
    bits_before = index & mask_before
    bits_after = index & mask_after
    
    # Shift bits after the position to make room for the new bit
    bits_after = bits_after << 1
    
    # Insert the new bit
    bit_mask = bit << position
    
    # Combine all parts
    return bits_before | bit_mask | bits_after
end

"""
    add_qubits(dms::AbstractMatrix{<:Complex}, num_qubits::Int) -> DensityMatrixSimulator

Add `num_qubits` qubits to the density matrix simulator, initializing them in the |0⟩ state.
Returns a new simulator with the expanded state.

# Arguments
- `dms`: The density matrix simulator to expand
- `num_qubits`: The number of qubits to add

# Returns
- A new density matrix with the expanded state
"""
function add_qubits(matrix::AbstractMatrix{<:Complex}, num_qubits::Int)
    # Validate input
    num_qubits >= 0 || throw(ArgumentError("Number of qubits to add must be non-negative"))
    
    # If no qubits to add, just return a copy
    if num_qubits == 0
        return copy(matrix)
    end
    
    # Handle empty density matrix case - initialize with |0⟩⟨0| state
    if isempty(matrix) || size(matrix, 1) == 0
        # Create a density matrix with 2^num_qubits x 2^num_qubits elements
        T = eltype(matrix)
        new_size = 2^num_qubits
        new_matrix = zeros(T, new_size, new_size)
        new_matrix[1, 1] = one(T)  # Set |0⟩⟨0| state
        return new_matrix
    end
    
    # Calculate new dimensions
    old_size = size(matrix, 1)
    new_size = old_size * (2^num_qubits)
    
    # Create new density matrix
    T = eltype(matrix)
    new_matrix = similar(matrix, new_size, new_size)
    fill!(new_matrix, zero(T))
    
    # Copy the old density matrix into the new one
    # Only the basis states where the new qubits are all |0⟩ get the original values
    for i in 0:(old_size-1)
        for j in 0:(old_size-1)
            new_matrix[i+1, j+1] = matrix[i+1, j+1]
        end
    end
    
    # Create new simulator with expanded state
    return new_matrix
end
