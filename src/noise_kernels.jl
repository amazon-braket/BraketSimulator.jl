@inline function get_amps_and_qubits(density_mat::AbstractDensityMatrix, qubits::Int...)
    n_amps   = size(density_mat, 1)
    n_qubits = Int(log2(n_amps))
    return n_amps, endian_qubits(n_qubits, qubits...)
end

function apply_noise!(n::PhaseFlip, dm::AbstractDensityMatrix{T}, qubit::Int) where {T}
    # K₀ = √(1.0 - n.probability) * I
    # K₁ = √(n.probability) * Z
    # ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
    n_amps, endian_qubit = get_amps_and_qubits(dm, qubit)
    p = n.probability
    Threads.@threads for idx = 0:length(dm)-1
        ix      = mod(idx, n_amps)
        jx      = div(idx, n_amps)
        i_bit   = ((1 << endian_qubit) & ix) >> endian_qubit
        j_bit   = ((1 << endian_qubit) & jx) >> endian_qubit
        i_value = (1.0 - 2.0 * i_bit)
        j_value = (1.0 - 2.0 * j_bit)
        factor  = (1.0 - p) + i_value * j_value * p
        @inbounds dm[CartesianIndex(ix + 1, jx + 1)] *= factor
    end
end

# K₀ = √(1.0 - n.probability) * I
# K₁ = √(n.probability / 3.0) * X
# K₂ = √(n.probability / 3.0) * Y
# K₃ = √(n.probability / 3.0) * Z
# ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
kraus_rep(n::Depolarizing) = (p = n.probability; return Kraus([
    √(1.0 - p) * matrix_rep(I()),
    √(p / 3.0) * matrix_rep(X()),
    √(p / 3.0) * matrix_rep(Y()),
    √(p / 3.0) * matrix_rep(Z()),
   ]))
kraus_rep(n::BitFlip) = (p = n.probability; return Kraus([√(1 - p) * matrix_rep(I()), √p * matrix_rep(X())]))
kraus_rep(n::PauliChannel) = Kraus([
    √(1 - n.probX - n.probY - n.probZ) * matrix_rep(I()),
    √n.probX * matrix_rep(X()),
    √n.probY * matrix_rep(Y()),
    √n.probZ * matrix_rep(Z()),
])
kraus_rep(n::AmplitudeDamping) = (γ = n.gamma; return Kraus([complex([1.0 0.0; 0.0 √(1.0 - γ)]), complex([0.0 √γ; 0.0 0.0])]))
kraus_rep(n::GeneralizedAmplitudeDamping) = (p = n.probability; γ = n.gamma; Kraus([
    √p * [1.0 0.0; 0.0 √(1.0 - γ)],
    √p * [0.0 √γ; 0.0 0.0],
    √(1.0 - p) * [√(1.0 - γ) 0.0; 0.0 1.0],
    √(1.0 - p) * [0.0 0.0; √γ 0.0],
   ]))
kraus_rep(n::PhaseDamping) = (γ = n.gamma; Kraus([complex([1.0 0.0; 0.0 √(1.0 - γ)]), complex([0.0 0.0; 0.0 √γ])]))

function kraus_rep(n::TwoQubitDepolarizing)
    I = diagm(ones(ComplexF64, 2))
    Z = diagm(ComplexF64[1.0; -1.0])
    X = ComplexF64[0.0 1.0; 1.0 0.0]
    Y = ComplexF64[0.0 -im; im 0.0]
    ks     = (I, X, Y, Z)
    p      = n.probability
    factor = √(p / 15.0)
    Ks     = [kron(ki, kj) for ki in ks, kj in ks]
    Ks[1] .*= √(1.0 - p)
    Ks[2:end] .*= factor
    return Kraus(vec(Ks))
end

function kraus_rep(n::TwoQubitPauliChannel)
    I = diagm(ones(ComplexF64, 2))
    Z = diagm(ComplexF64[1.0; -1.0])
    X = ComplexF64[0.0 1.0; 1.0 0.0]
    Y = ComplexF64[0.0 -im; im 0.0]
    k_dict = Dict('I' => I, 'X' => X, 'Y' => Y, 'Z' => Z)
    total_p = sum(values(n.probabilities))
    Ks = [diagm(√(1.0 - total_p) * ones(ComplexF64, 4))]
    for (k, v) in n.probabilities
        k != "II" && push!(Ks, √v * kron(k_dict[k[1]], k_dict[k[2]]))
    end
    return Kraus(Ks)
end

function kraus_rep(n::TwoQubitDephasing)
    # K₀ = √(1.0 - n.probability) * II 
    # K₁ = √(n.probability / 3.0) * IZ
    # K₂ = √(n.probability / 3.0) * ZI
    # K₃ = √(n.probability / 3.0) * ZZ
    # ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
    I  = diagm(ones(ComplexF64, 2))
    Z  = diagm(ComplexF64[1.0; -1.0])
    II = kron(I, I)
    IZ = kron(I, Z)
    ZI = kron(Z, I)
    ZZ = kron(Z, Z)
    p  = n.probability
    Ks = [
        √(1-p)   * II,
        √(p/3.0) * IZ,
        √(p/3.0) * ZI,
        √(p/3.0) * ZZ,
    ]
    return Kraus(Ks)
end

apply_noise!(n::N, dm::S, qubits::Int...) where {T,S<:AbstractDensityMatrix{T},N<:Noise} =
    apply_noise!(kraus_rep(n), dm, qubits...)

function apply_noise!(n::Kraus, dm::AbstractDensityMatrix{T}, qubit::Int) where {T}
    k_mats = ntuple(ix -> SMatrix{2,2,ComplexF64}(n.matrices[ix]), length(n.matrices))
    k_mats_conj =
        ntuple(ix -> SMatrix{2,2,ComplexF64}(adjoint(n.matrices[ix])), length(n.matrices))
    # ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
    n_amps, endian_qubit = get_amps_and_qubits(dm, qubit)
    Threads.@threads for ix = 0:div(n_amps, 2)-1
        # maybe not diagonal
        lower_ix  = pad_bit(ix, endian_qubit)
        higher_ix = flip_bit(lower_ix, endian_qubit) + 1
        lower_ix += 1
        @inbounds begin
            ρ_00 = dm[CartesianIndex(lower_ix, lower_ix)]
            ρ_01 = dm[CartesianIndex(lower_ix, higher_ix)]
            ρ_10 = dm[CartesianIndex(higher_ix, lower_ix)]
            ρ_11 = dm[CartesianIndex(higher_ix, higher_ix)]
            sm_ρ = SMatrix{2,2,ComplexF64}(ρ_00, ρ_10, ρ_01, ρ_11)

            k_ρ = k_mats[1] * sm_ρ * k_mats_conj[1]
            for mat_ix = 2:length(k_mats)
                k_ρ += k_mats[mat_ix] * sm_ρ * k_mats_conj[mat_ix]
            end
            dm[CartesianIndex(lower_ix, lower_ix)]   = k_ρ[1, 1]
            dm[CartesianIndex(lower_ix, higher_ix)]  = k_ρ[1, 2]
            dm[CartesianIndex(higher_ix, lower_ix)]  = k_ρ[2, 1]
            dm[CartesianIndex(higher_ix, higher_ix)] = k_ρ[2, 2]
        end
    end
end

function apply_noise!(n::Kraus, dm::AbstractDensityMatrix{T}, target1::Int, target2::Int) where {T}
    k_mats = ntuple(ix -> SMatrix{4,4,ComplexF64}(n.matrices[ix]), length(n.matrices))
    k_mats_conj =
        ntuple(ix -> SMatrix{4,4,ComplexF64}(adjoint(n.matrices[ix])), length(n.matrices))
    # ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
    n_amps, (endian_t1, endian_t2) = get_amps_and_qubits(dm, target1, target2)
    small_t, big_t = minmax(endian_t1, endian_t2)
    Threads.@threads for ix = 0:div(n_amps, 4)-1
        # maybe not diagonal
        padded_ix = pad_bits(ix, (small_t, big_t))
        ix_00 = padded_ix + 1
        ix_01 = flip_bit(padded_ix, endian_t1) + 1
        ix_10 = flip_bit(padded_ix, endian_t2) + 1
        ix_11 = flip_bits(padded_ix, (endian_t1, endian_t2)) + 1
        @inbounds begin
            ixs =
                CartesianIndex.(
                    collect(
                        Iterators.product(
                            (ix_00, ix_10, ix_01, ix_11),
                            (ix_00, ix_10, ix_01, ix_11),
                        ),
                    )
                )
            dm_vec = view(dm, ixs)
            ρ   = SMatrix{4,4,ComplexF64}(dm_vec)
            k_ρ = k_mats[1] * ρ * k_mats_conj[1]
            for mat_ix = 2:length(k_mats)
                k_ρ += k_mats[mat_ix] * ρ * k_mats_conj[mat_ix]
            end
            dm_vec[:] = k_ρ[:]
        end
    end
end

function apply_noise!(n::Kraus, dm::AbstractDensityMatrix{T}, targets::Int...) where {T}
    k_mats = n.matrices
    k_mats_conj = ntuple(ix -> adjoint(n.matrices[ix]), length(n.matrices))
    # ρ = ∑ᵢ Kᵢ ρ Kᵢ\^†
    n_amps, endian_ts = get_amps_and_qubits(dm, targets...)
    ordered_ts = sort(collect(endian_ts))
    n_targets  = length(targets)
    flip_list  = map(0:2^n_targets-1) do target_amp
        f_vals = Bool[(((1 << target) & target_amp) >> target) for target = 0:n_targets-1]
        return ordered_ts[f_vals]
    end
    Threads.@threads for ix = 0:div(n_amps, 2^n_targets)-1
        padded_ix = pad_bits(ix, ordered_ts)
        ixs = map(flip_list) do to_flip
            flip_bits(padded_ix, to_flip) + 1
        end
        ix_pairs = CartesianIndex.(collect(Iterators.product(ixs, ixs)))
        @views @inbounds begin
            ρ   = dm[ix_pairs]
            k_ρ = k_mats[1] * ρ * k_mats_conj[1]
            for mat_ix = 2:length(k_mats)
                k_ρ += k_mats[mat_ix] * ρ * k_mats_conj[mat_ix]
            end
            dm[ix_pairs] = k_ρ[:, :]
        end
    end
end
