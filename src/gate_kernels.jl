@inline function pad_bit(amp_index::Ti, t::Tj)::Ti where {Ti<:Integer,Tj<:Integer}
    left = (amp_index >> t) << t
    right = amp_index - left
    return (left << one(Ti)) ⊻ right
end
function pad_bits(ix::Ti, to_pad)::Ti where {Ti<:Integer}
    padded_ix = ix
    for bit in to_pad
        padded_ix = pad_bit(padded_ix, bit)
    end
    return padded_ix
end

@inline function flip_bit(amp_index::Ti, t::Tj)::Ti where {Ti<:Integer,Tj<:Integer}
    return amp_index ⊻ (one(Ti) << t)
end
function flip_bits(ix::Ti, to_flip)::Ti where {Ti<:Integer}
    flipped_ix = ix
    for f_val in to_flip
        flipped_ix = flip_bit(flipped_ix, f_val)
    end
    return flipped_ix
end
@inline endian_qubits(nq::Int, qubit::Int) = nq - 1 - qubit
@inline endian_qubits(nq::Int, qubits::Int...) = nq .- 1 .- qubits
@inline function get_amps_and_qubits(state_vec::AbstractStateVector, qubits::Int...)
    n_amps = length(state_vec)
    nq = Int(log2(n_amps))
    return n_amps, endian_qubits(nq, qubits...)
end

matrix_rep(g::H) = SMatrix{2,2}(complex((1 / √2) * [1.0 1.0; 1.0 -1.0]))
matrix_rep(g::X) = SMatrix{2,2}(complex([0.0 1.0; 1.0 0.0]))
matrix_rep(g::Y) = SMatrix{2,2}(complex([0.0 -im; im 0.0]))
matrix_rep(g::Z) = SMatrix{2,2}(complex([1.0 0.0; 0.0 -1.0]))
matrix_rep(g::I) = SMatrix{2,2}(complex([1.0 0.0; 0.0 1.0]))
matrix_rep(g::V) = SMatrix{2,2}(0.5 * [1.0+im 1.0-im; 1.0-im 1.0+im])
matrix_rep(g::Vi) = SMatrix{2,2}(0.5 * [1.0-im 1.0+im; 1.0+im 1.0-im])
matrix_rep(g::S) = SMatrix{2,2}([1.0 0.0; 0.0 im])
matrix_rep(g::Si) = SMatrix{2,2}([1.0 0.0; 0.0 -im])
matrix_rep(g::T) = SMatrix{2,2}([1.0 0.0; 0.0 exp(im * π / 4.0)])
matrix_rep(g::Ti) = SMatrix{2,2}([1.0 0.0; 0.0 exp(-im * π / 4.0)])

matrix_rep(g::CNot) = SMatrix{4,4}(complex([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0]))
matrix_rep(g::CY) = SMatrix{4,4}([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 -im; 0.0 0.0 im 0.0])
matrix_rep(g::CZ) = SMatrix{4,4}(complex([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0]))
matrix_rep(g::CV) = SMatrix{4,4}([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.5+0.5im 0.5-0.5im; 0.0 0.0 0.5-0.5im 0.5+0.5im])
matrix_rep(g::CCNot) = SMatrix{8,8}(complex([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0]))
matrix_rep(g::CSwap) = SMatrix{8,8}(complex([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]))
matrix_rep(g::PRx) = SMatrix{2,2}(
    [
         cos(g.angle[1] / 2.0) -im*exp(-im*g.angle[2])*sin(g.angle[1]/2.0)
        -im*exp(im*g.angle[2])*sin(g.angle[1]/2.0) cos(g.angle[1] / 2.0)
    ],
)
matrix_rep(g::Rz) =
    SMatrix{2,2}([exp(-im * g.angle[1] / 2.0) 0.0; 0.0 exp(im * g.angle[1] / 2.0)])
matrix_rep(g::Rx) = SMatrix{2,2}(
    [
        cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0)
        -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0)
    ],
)
matrix_rep(g::Ry) = SMatrix{2,2}(
    complex(
        [
            cos(g.angle[1] / 2.0) -sin(g.angle[1] / 2.0)
            sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0)
        ],
    ),
)
matrix_rep(g::GPi) =
    SMatrix{2,2}(complex([0 exp(-im * g.angle[1]); exp(im * g.angle[1]) 0]))
matrix_rep(g::GPi2) =
    SMatrix{2,2}(1/√2 * complex([1.0 -im*exp(-im * g.angle[1]); -im*exp(im * g.angle[1]) 1.0]))
matrix_rep(g::MS) = SMatrix{4,4}(
    complex(
        [
            cos(g.angle[3] / 2.0) 0.0 0.0 -im*exp(-im * (g.angle[1] + g.angle[2]))*sin(g.angle[3] / 2)
            0.0 cos(g.angle[3] / 2) -im*exp(-im * (g.angle[1] - g.angle[2]))*sin(g.angle[3] / 2) 0.0
            0.0 -im*exp(im * (g.angle[1] - g.angle[2]))*sin(g.angle[3] / 2) cos(g.angle[3] / 2) 0.0
            -im*exp(im * (g.angle[1] + g.angle[2]))*sin(g.angle[3] / 2) 0.0 0.0 cos(g.angle[3] / 2)
        ],
    ),
)
matrix_rep(g::PhaseShift) = SMatrix{2,2}([1.0 0.0; 0.0 exp(im * g.angle[1])])
matrix_rep(g::XX) = SMatrix{4,4}(
    [
        cos(g.angle[1] / 2.0) 0.0 0.0 -im*sin(g.angle[1] / 2)
        0.0 cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0) 0.0
        0.0 -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        -im*sin(g.angle[1] / 2.0) 0.0 0.0 cos(g.angle[1] / 2.0)
    ],
)
matrix_rep(g::XY) = SMatrix{4,4}(
    [
        1.0 0.0 0.0 0.0
        0.0 cos(g.angle[1] / 2.0) im*sin(g.angle[1] / 2.0) 0.0
        0.0 im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        0.0 0.0 0.0 1.0
    ],
)
matrix_rep(g::YY) = SMatrix{4,4}(
    [
        cos(g.angle[1] / 2.0) 0.0 0.0 im*sin(g.angle[1] / 2)
        0.0 cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0) 0.0
        0.0 -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        im*sin(g.angle[1] / 2.0) 0.0 0.0 cos(g.angle[1] / 2.0)
    ],
)
matrix_rep(g::ZZ) = SMatrix{4,4}(
    diagm([
        exp(-im * g.angle[1] / 2.0),
        exp(im * g.angle[1] / 2.0),
        exp(im * g.angle[1] / 2.0),
        exp(-im * g.angle[1] / 2.0),
    ]),
)
# 1/√2 * (IX - XY)
matrix_rep(g::ECR) = SMatrix{4,4}(1/√2 * [0 1 0 im; 1 0 -im 0; 0 im 0 1; -im 0 1 0])
matrix_rep(g::Unitary) = g.matrix

apply_gate!(::Val{false}, g::I, state_vec::StateVector{T}, qubits::Int...) where {T<:Complex} =
    return
apply_gate!(::Val{true}, g::I, state_vec::StateVector{T}, qubits::Int...) where {T<:Complex} =
    return

function apply_gate!(
    g_mat::Union{SMatrix{2,2,T}, Diagonal{T,SVector{2,T}}},
    state_vec::StateVector{T},
    qubit::Int,
) where {T<:Complex}
    n_amps, endian_qubit = get_amps_and_qubits(state_vec, qubit)
    n_tasks = n_amps >> 1
    n_chunks = max(div(n_tasks, CHUNK_SIZE), 1)
    flipper = 1 << endian_qubit
    is_small_target = flipper < CHUNK_SIZE
    g_00, g_10, g_01, g_11 = g_mat
    Threads.@threads for c_ix = 0:n_chunks-1
        my_amps = if n_chunks > 1
            c_ix*CHUNK_SIZE:((c_ix+1)*CHUNK_SIZE-1)
        else
            0:n_tasks-1
        end
        lower_ix = pad_bit(my_amps[1], endian_qubit) + 1
        higher_ix = lower_ix + flipper
        for task_amp = 0:length(my_amps)-1
            if is_small_target && div(task_amp, flipper) > 0 && mod(task_amp, flipper) == 0
                lower_ix = higher_ix
                higher_ix = lower_ix + flipper
            end
            @inbounds begin
                lower_amp = state_vec[lower_ix]
                higher_amp = state_vec[higher_ix]
                state_vec[lower_ix] = g_00 * lower_amp + g_01 * higher_amp
                state_vec[higher_ix] = g_10 * lower_amp + g_11 * higher_amp
            end
            lower_ix += 1
            higher_ix += 1
        end
    end
end

# generic two-qubit non controlled unitaries
function apply_gate!(
    g_mat::M,
    state_vec::StateVector{T},
    t1::Int,
    t2::Int,
) where {T<:Complex,M<:Union{SMatrix{4,4,T},Diagonal{T,SVector{4,T}}}}
    n_amps, (endian_t1, endian_t2) = get_amps_and_qubits(state_vec, t1, t2)
    small_t, big_t = minmax(endian_t1, endian_t2)
    chunk_size     = CHUNK_SIZE
    n_tasks        = n_amps >> 2
    n_chunks       = max(div(n_tasks, CHUNK_SIZE), 1)
    chunked_amps   = collect(Iterators.partition(0:n_tasks-1, CHUNK_SIZE))
    Threads.@threads for c_ix = 1:n_chunks
        for ix in chunked_amps[c_ix]
            ix_00   = pad_bit(pad_bit(ix, small_t), big_t)
            ix_10   = flip_bit(ix_00, endian_t2)
            ix_01   = flip_bit(ix_00, endian_t1)
            ix_11   = flip_bit(ix_10, endian_t1)
            ind_vec = SVector(ix_00 + 1, ix_01 + 1, ix_10 + 1, ix_11 + 1)
            @inbounds begin
                amps = SVector{4,T}(
                    state_vec[ix_00+1],
                    state_vec[ix_01+1],
                    state_vec[ix_10+1],
                    state_vec[ix_11+1],
                )
                new_amps = g_mat * amps
                state_vec[ind_vec[1]] = new_amps[1]
                state_vec[ind_vec[2]] = new_amps[2]
                state_vec[ind_vec[3]] = new_amps[3]
                state_vec[ind_vec[4]] = new_amps[4]
            end
        end
    end
end

for (V, f) in ((false, :identity), (true, :conj))
    @eval begin
        apply_gate!(::Val{$V},
                    g::G,
                    state_vec::StateVector{T},
                    qubits::Int...) where {G<:Gate,T<:Complex} = apply_gate!($f(matrix_rep(g)), state_vec, qubits...)
        apply_gate!(::Val{$V},
                    g::Unitary,
                    state_vec::StateVector{T},
                    qubit::Int) where {T<:Complex} = apply_gate!(SMatrix{2,2,T}($f(matrix_rep(g))), state_vec, qubit)
        apply_gate!(::Val{$V},
                    g::Unitary,
                    state_vec::StateVector{T},
                    q1::Int,
                    q2::Int) where {T<:Complex} = apply_gate!(SMatrix{4,4,T}($f(matrix_rep(g))), state_vec, q2, q1)
    end
end

apply_gate!(g::G, args...) where {G<:Gate} = apply_gate!(Val(false), g, args...)

# controlled unitaries
for (cph, ind) in
    ((:CPhaseShift, 4), (:CPhaseShift00, 1), (:CPhaseShift01, 3), (:CPhaseShift10, 2))
    @eval begin
        matrix_rep(g::$cph) = Diagonal(
            SVector{4,ComplexF64}(
                setindex!(ones(ComplexF64, 4), exp(im * g.angle[1]), $ind),
            ),
        )
    end
end

for (sw, factor) in ((:Swap, 1.0), (:ISwap, im), (:PSwap, :(exp(im * g.angle[1]))))
    @eval begin
        matrix_rep(g::$sw) = SMatrix{4,4,ComplexF64}(
            [1.0 0.0 0.0 0.0 0.0 0.0 $factor 0.0 0.0 $factor 0.0 0.0 0.0 0.0 0.0 1.0],
        )
    end
end

function apply_controlled_gate!(
    g_matrix::SMatrix{2, 2, ComplexF64},
    c_bit::Bool,
    state_vec::StateVector{T},
    control::Int,
    target::Int,
) where {T<:Complex}
    n_amps, (endian_control, endian_target) =
        get_amps_and_qubits(state_vec, control, target)
    
    small_t, big_t = minmax(endian_control, endian_target)
    g_00, g_10, g_01, g_11 = g_matrix
    Threads.@threads for ix = 0:div(n_amps, 4)-1
        ix_00 = pad_bit(pad_bit(ix, small_t), big_t)
        ix_10 = flip_bit(ix_00, endian_control)
        ix_01 = flip_bit(ix_00, endian_target)
        ix_11 = flip_bit(ix_01, endian_control)
        lower_ix  = c_bit ? ix_10 + 1 : ix_00 + 1
        higher_ix = c_bit ? ix_11 + 1 : ix_01 + 1
        @inbounds begin
            lower_amp  = state_vec[lower_ix]
            higher_amp = state_vec[higher_ix]
            state_vec[lower_ix]  = g_00 * lower_amp + g_01 * higher_amp
            state_vec[higher_ix] = g_10 * lower_amp + g_11 * higher_amp
        end
    end
    return
end
function apply_controlled_gate!(
    g_matrix::SMatrix{4, 4, ComplexF64},
    c_bit::Bool,
    state_vec::StateVector{T},
    control::Int,
    t1::Int,
    t2::Int,
) where {T<:Complex}
    n_amps, (endian_control, endian_t1, endian_t2) =
        get_amps_and_qubits(state_vec, control, t1, t2)
    small_t, mid_t, big_t = sort([endian_control, endian_t1, endian_t2])
    Threads.@threads for ix = 0:div(n_amps, 8)-1
        ix_00 = pad_bits(ix, (small_t, mid_t, big_t))
        if c_bit
            ix_00 = flip_bit(ix_00, endian_control)
        end
        ix_10 = flip_bit(ix_00, endian_t2)
        ix_01 = flip_bit(ix_00, endian_t1)
        ix_11 = flip_bit(ix_01, endian_t2)
        ix_vec = SVector{4,Int}(ix_00 + 1, ix_01 + 1, ix_10 + 1, ix_11 + 1)
        @views @inbounds begin
            amps = SVector{4,T}(state_vec[ix_vec])
            state_vec[ix_vec] = g_matrix * amps
        end
    end
    return
end
# doubly controlled unitaries
function apply_controlled_gate!(
    g_matrix::SMatrix{2, 2, ComplexF64},
    c1_bit::Bool,
    c2_bit::Bool,
    state_vec::StateVector{T},
    c1::Int,
    c2::Int,
    t::Int,
) where {T<:Complex}
    n_amps, (endian_c1, endian_c2, endian_target) =
        get_amps_and_qubits(state_vec, c1, c2, t)
    small_t, mid_t, big_t = sort([endian_c1, endian_c2, endian_target])
    g_00, g_10, g_01, g_11 = g_matrix
    Threads.@threads for ix = 0:div(n_amps, 8)-1
        # insert 0 at c1, 0 at c2, 0 at target
        padded_ix = pad_bits(ix, [small_t, mid_t, big_t])
        # flip c1 and c2 if needed
        lower_ix = padded_ix
        for (c_val, c_bit) in zip((c1_bit, c2_bit), (endian_c1, endian_c2))
            c_val == 1 && (lower_ix = flip_bit(lower_ix, c_bit))
        end
        # flip target
        higher_ix = flip_bit(lower_ix, endian_target) + 1
        lower_ix += 1
        @inbounds begin
            lower_amp  = state_vec[lower_ix]
            higher_amp = state_vec[higher_ix]
            state_vec[lower_ix]  = g_00 * lower_amp + g_01 * higher_amp
            state_vec[higher_ix] = g_10 * lower_amp + g_11 * higher_amp
        end
    end
    return
end
for (V, f) in ((true, :conj), (false, :identity))
    @eval begin
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{1},
            g::G,
            tg::TG,
            state_vec::StateVector{T},
            control_bits::NTuple{1,Int},
            control::Int,
            target::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(tg)), control_bits[1] == 1, state_vec, control, target)
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{1},
            g::G,
            tg::TG,
            state_vec::StateVector{T},
            control_bits::NTuple{1,Int},
            control::Int,
            t1::Int,
            t2::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(tg)), control_bits[1] == 1, state_vec, control, t1, t2)
        # doubly controlled unitaries
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{2},
            g::G,
            tg::TG,
            state_vec::StateVector{T},
            control_bits::NTuple{2,Int},
            c1::Int,
            c2::Int,
            t::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(tg)), control_bits[1] == 1, control_bits[2] == 1, state_vec, c1, c2, t)
    end
end

for (cg, tg, nc) in (
        (:CNot, :X, 1),
        (:CY, :Y, 1),
        (:CZ, :Z, 1),
        (:CV, :V, 1),
        (:CSwap, :Swap, 1),
        (:CCNot, :X, 2),
    ),
    (Vc, f) in ((false, :identity), (true, :conj))

    @eval begin
        apply_gate!(
            ::Val{$Vc},
            g::$cg,
            state_vec::StateVector{T},
            qubits::Int...,
        ) where {T<:Complex} =
            apply_controlled_gate!(Val($Vc), Val($nc), g, $tg(), state_vec, ntuple(i->1, Val($nc)), qubits...)
    end
end


function apply_gate!(
    factor::Complex,
    g::Control{MultiQubitPhaseShift{N}, B},
    state_vec::StateVector{T},
    qubits::Int...,
) where {T<:Complex, N, B}
    bvs   = g.bitvals
    g_mat = ones(ComplexF64, 2^N)
    ix    = 0
    for (bi, bv) in enumerate(bvs)
        ix += bv << (bi - 1)
    end
    ix += 1
    g_mat[ix] = exp(factor*g.g.angle[1])
    apply_gate!(Diagonal(SVector{2^N, ComplexF64}(g_mat)), state_vec, qubits...)
end

# arbitrary number of targets
for (V, f) in ((false, :identity), (true, :conj))
    @eval begin
        apply_gate!(
            ::Val{$V},
            g::Control{G, B},
            state_vec::StateVector{T},
            qubits::Int...,
        ) where {T<:Complex, G<:Gate, B} =
            apply_controlled_gate!(Val($V), Val(B), g, g.g, state_vec, g.bitvals, qubits...)
        apply_gate!(
            ::Val{$V},
            g::Control{MultiQubitPhaseShift{N}, B},
            state_vec::StateVector{T},
            qubits::Int...,
        ) where {T<:Complex, N, B} = apply_gate!($f(im), g, state_vec, qubits...)
        apply_gate!(
            ::Val{$V},
            g::Unitary,
            state_vec::StateVector{T},
            ts::Vararg{Int,NQ},
        ) where {T<:Complex,NQ} = 
            apply_gate!($f(SMatrix{2^NQ, 2^NQ, ComplexF64}(matrix_rep(g))), state_vec, ts...)
    end
end
function apply_gate!(
    g_mat::Union{SMatrix{N, N, T}, Diagonal{T, SVector{N, T}}},
    state_vec::StateVector{T},
    ts::Vararg{Int,NQ},
) where {T<:Complex,NQ,N}
    n_amps, endian_ts = get_amps_and_qubits(state_vec, ts...)
    endian_ts  = tuple(endian_ts...)
    ordered_ts = sort(collect(endian_ts))
    flip_list  = map(0:N-1) do t
        f_vals = Bool[(((1 << f_ix) & t) >> f_ix) for f_ix = 0:NQ-1]
        return ordered_ts[f_vals]
    end
    Threads.@threads for ix = 0:div(n_amps, N)-1
        padded_ix = pad_bits(ix, ordered_ts)
        ixs = SVector{N,Int}(flip_bits(padded_ix, f) + 1 for f in flip_list)
        @views @inbounds begin
            amps = SVector{N,T}(state_vec[ixs])
            new_amps = g_mat * amps
            state_vec[ixs] = new_amps
        end
    end
    return
end