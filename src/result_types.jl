# In-place version of alias table building algorithm from 
# `StatsBase.jl`. We use this so that we can re-use the
# arguments `a`, `alias`, `larges`, and `smalls` **without**
# reallocating them, which is useful if we are using the same
# simulator instantiation for many circuit simulations (e.g.
# in a large batch). Based on the [alias method](https://en.wikipedia.org/wiki/Alias_method)
# for sampling from discrete probability distributions.
function make_alias_table!(
    weights::AbstractVector,
    wsum,
    acceptance_probs::AbstractVector{Float64},
    alias::AbstractVector{Int},
    larges::AbstractVector{Int},
    smalls::AbstractVector{Int},
)
    n = length(weights)
    length(acceptance_probs) == length(alias) == n ||
        throw(DimensionMismatch("Inconsistent array lengths. length(acceptance_probs) = $(length(acceptance_probs)), length(alias) = $(length(alias)), n = $n"))

    ac = n / wsum
    acceptance_probs .= weights .* ac

    kl = 0  # actual number of larges
    ks = 0  # actual number of smalls

    @inbounds for i = 1:n
        ai = acceptance_probs[i]
        if ai > 1.0
            larges[kl+=1] = i  # push to larges
        elseif ai < 1.0
            smalls[ks+=1] = i  # push to smalls
        end
    end

    @inbounds while kl > 0 && ks > 0
        s = smalls[ks]
        ks -= 1  # pop from smalls
        l = larges[kl]
        kl -= 1  # pop from larges
        alias[s] = l
        al = acceptance_probs[l] = (acceptance_probs[l] - 1.0) + acceptance_probs[s]
        if al > 1.0
            larges[kl+=1] = l  # push to larges
        else
            smalls[ks+=1] = l  # push to smalls
        end
    end

    # this loop should be redundant, except for rounding
    for i = 1:ks
        @inbounds acceptance_probs[smalls[i]] = 1.0
    end
    nothing
end

function samples(simulator::AbstractSimulator)
    sim_probabilities = probabilities(simulator)
    weight_vector = Weights(sim_probabilities)
    n_amplitudes  = 2^simulator.qubit_count
    rng = Random.Xoshiro()
    if simulator.qubit_count < 30 # build alias tables etc
        inds   = 0:(n_amplitudes-1)
        ap     = simulator._ap
        alias  = simulator._alias
        larges = simulator._larges
        smalls = simulator._smalls
        make_alias_table!(weight_vector, sum(weight_vector), ap, alias, larges, smalls)
        sampler = Random.Sampler(rng, 1:n_amplitudes)
        for shot_index = 1:simulator.shots
            amplitude_index = rand(rng, sampler)
            simulator.shot_buffer[shot_index] = rand(rng) < ap[amplitude_index] ? amplitude_index - 1 : alias[amplitude_index] - 1
        end
    else
        # Direct sampling
        for shot_index = 1:simulator.shots
            simulator.shot_buffer[shot_index] = StatsBase.sample(rng, weight_vector) - 1
        end
    end
    return simulator.shot_buffer
end

calculate(sv::Braket.StateVector, sim::AbstractSimulator) = state_vector(sim)
function calculate(amplitude::Braket.Amplitude, sim::AbstractSimulator)
    state = collect(state_vector(sim))
    rev_states = reverse.(amplitude.states)
    state_ints = [
        sum(tryparse(Int, string(amp_state[qubit])) * 2^(qubit - 1) for qubit = 1:length(amp_state)) for amp_state in rev_states
    ]
    return Dict(
        basis_state => state[state_int+1] for
        (basis_state, state_int) in zip(amplitude.states, state_ints)
    )
end

# need to handle reordered targets here
function marginal_probability(probs::Vector{T}, qubit_count::Int, targets) where {T<:Real}
    unused_qubits = setdiff(collect(0:qubit_count-1), targets)
    endian_unused = qubit_count .- unused_qubits .- 1
    final_probs   = zeros(Float64, 2^length(targets))
    qubit_combos  = vcat([Int[]], collect(combinations(endian_unused)))
    Threads.@threads for ix = 0:2^length(targets)-1
        padded_ix = ix
        for pad_q in sort(endian_unused)
            padded_ix = pad_bit(padded_ix, pad_q)
        end
        # here we generate for the **output** ix all amplitude indices in the
        # **full** 1:2^n_qubits list that will go into the marginal sum
        flipped_inds = Vector{Int64}(undef, length(qubit_combos))
        for (c_ix, flipped_qubits) in enumerate(qubit_combos)
            flipped_ix = padded_ix
            for flip_qubit in flipped_qubits
                flipped_ix = flip_bit(flipped_ix, flip_qubit)
            end
            flipped_inds[c_ix] = flipped_ix + 1
        end
        @views begin
            @inbounds sum_val = sum(probs[flipped_inds])
            final_probs[ix+1] = sum_val
        end
    end
    return final_probs
end

function permute_probability(probabilities::Vector{T}, qubit_count::Int, targets) where {T<:Real}
    new_probabilities = zeros(T, 2^length(targets))
    Threads.@threads for ix = 0:2^length(targets)-1
        unpermed_bits = [(ix >> qubit) & 1 for qubit in 0:qubit_count - 1]
        permed_bits   = unpermed_bits[qubit_count .- targets] 
        new_ix = 0
        for (permuted_bit, qubit) in zip(permed_bits, 0:qubit_count - 1)
            new_ix += permuted_bit << qubit
        end
        new_probabilities[new_ix+1] = probabilities[ix+1]
    end
    return new_probabilities
end

function permute_density_matrix(ρ::Matrix{T}, qubit_count::Int, targets) where {T}
    new_ρ = zeros(T, 2^length(targets), 2^length(targets))
    Threads.@threads for ix = 0:2^length(targets)-1
        unpermed_ix_bits = [(ix >> qubit) & 1 for qubit in 0:qubit_count - 1]
        permed_ix_bits   = unpermed_ix_bits[qubit_count .- targets] 
        new_ix = 0
        for (permuted_bit, qubit) in zip(permed_ix_bits, 0:qubit_count - 1)
            new_ix += permuted_bit << qubit
        end
        for jx = 0:2^length(targets)-1
            unpermed_jx_bits = [(jx >> qubit) & 1 for qubit in 0:qubit_count - 1]
            permed_jx_bits   = unpermed_jx_bits[qubit_count .- targets] 
            new_jx = 0
            for (permuted_bit, qubit) in zip(permed_jx_bits, 0:qubit_count - 1)
                new_jx += permuted_bit << qubit
            end
            new_ρ[new_ix+1, new_jx+1] = ρ[ix+1, jx+1]
        end
    end
    return new_ρ
end

function calculate(probability::Braket.Probability, sim::AbstractSimulator)
    targets  = probability.targets
    probs    = probabilities(sim)
    n_qubits = qubit_count(sim)
    (isempty(targets) || [targets...] == collect(0:n_qubits-1)) && return probs
    # reordered
    length(targets) == n_qubits && return permute_probability(probs, n_qubits, [targets...])
    return marginal_probability(probs, n_qubits, targets)
end

function calculate(expectation_result::Braket.Expectation, sim::AbstractSimulator)
    obs             = expectation_result.observable
    targets         = isempty(expectation_result.targets) ? collect(0:qubit_count(sim)-1) : expectation_result.targets
    obs_qubit_count = qubit_count(obs)
    length(targets) == obs_qubit_count && return expectation(sim, obs, targets...)
    return [expectation(sim, obs, target) for target in targets]
end
expectation_op_squared(sim, obs::Braket.Observables.StandardObservable, target::Int) = 1.0
expectation_op_squared(sim, obs::Braket.Observables.I, target::Int) = 1.0
function expectation_op_squared(sim, obs::Braket.Observables.TensorProduct, targets::Int...)
    all(
        factor isa Braket.Observables.StandardObservable || factor isa Braket.Observables.I for
        factor in obs.factors
    ) && return 1.0
    sq_factors = map(obs.factors) do factor
        (factor isa Braket.Observables.StandardObservable || factor isa Braket.Observables.I) &&
            return Braket.Observables.I()
        factor isa Braket.Observables.HermitianObservable &&
            return Braket.Observables.HermitianObservable(factor.matrix * factor.matrix)
    end
    sq_tensor_prod = Braket.Observables.TensorProduct(sq_factors)
    return expectation(sim, sq_tensor_prod, targets...)
end
function expectation_op_squared(
    sim,
    obs::Braket.Observables.HermitianObservable,
    targets::Int...,
)
    return expectation(
        sim,
        Braket.Observables.HermitianObservable(obs.matrix * obs.matrix),
        targets...,
    )
end

function apply_observable!(
    observable::Braket.Observables.TensorProduct,
    sv_or_dm::T,
    targets::Int...,
) where {T<:AbstractVecOrMat{<:Complex}}
    target_ix = 1
    for factor in observable.factors
        factor_n_qubits = qubit_count(factor)
        factor_targets =
            factor_n_qubits == 1 ? targets[target_ix] : targets[target_ix:target_ix+factor_n_qubits-1]
        target_ix += factor_n_qubits
        sv_or_dm = apply_observable!(factor, sv_or_dm, factor_targets...)
    end
    return sv_or_dm
end
apply_observable(
    observable::O,
    sv_or_dm,
    target::Int...,
) where {O<:Braket.Observables.Observable} =
    apply_observable!(observable, deepcopy(sv_or_dm), target...)

function calculate(variance::Braket.Variance, sim::AbstractSimulator)
    obs     = variance.observable
    targets = isnothing(variance.targets) ? collect(0:qubit_count(sim)-1) : variance.targets
    obs_qubit_count = qubit_count(obs)
    if length(targets) == obs_qubit_count
        var2 = expectation_op_squared(sim, obs, targets...)
        mean = expectation(sim, obs, targets...)
        return var2 - mean^2
    else
        return map(collect(targets)) do target
            var2 = expectation_op_squared(sim, obs, target)
            mean = expectation(sim, obs, target)
            return var2 - mean^2
        end
    end
end

function calculate(dm::Braket.DensityMatrix, sim::AbstractSimulator)
    ρ = density_matrix(sim)
    full_qubits = collect(0:qubit_count(sim)-1)
    (collect(dm.targets) == full_qubits || isempty(dm.targets)) && return ρ
    length(dm.targets) == sim.qubit_count && return permute_density_matrix(ρ, sim.qubit_count, collect(dm.targets))
    # otherwise must compute a partial trace
    return partial_trace(ρ, dm.targets)
end
