@inline function pad_bit(amp_index::Ti, bit::Tj)::Ti where {Ti<:Integer,Tj<:Integer}
    left  = (amp_index >> bit) << bit
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

@inline function flip_bit(amp_index::Ti, bit::Tj)::Ti where {Ti<:Integer,Tj<:Integer}
    return amp_index ⊻ (one(Ti) << bit)
end
function flip_bits(ix::Ti, to_flip)::Ti where {Ti<:Integer}
    flipped_ix = ix
    for bit in to_flip
        flipped_ix = flip_bit(flipped_ix, bit)
    end
    return flipped_ix
end
@inline endian_qubits(n_qubits::Int, qubit::Int) = n_qubits - 1 - qubit
@inline endian_qubits(n_qubits::Int, qubits::Int...) = n_qubits .- 1 .- qubits
@inline function get_amps_and_qubits(state_vec::AbstractStateVector, qubits::Int...)
    n_amps   = length(state_vec)
    n_qubits = Int(log2(n_amps))
    return n_amps, endian_qubits(n_qubits, qubits...)
end

matrix_rep_raw(::H) = SMatrix{2,2}(complex((1 / √2) * [1.0 1.0; 1.0 -1.0]))
matrix_rep_raw(::X) = SMatrix{2,2}(complex([0.0 1.0; 1.0 0.0]))
matrix_rep_raw(::Y) = SMatrix{2,2}(complex([0.0 -im; im 0.0]))
matrix_rep_raw(::Z) = SMatrix{2,2}(complex([1.0 0.0; 0.0 -1.0]))
matrix_rep_raw(g::I, qc::Int=1) = SMatrix{2^qc,2^qc}(complex(diagm(ones(2^qc))))
matrix_rep_raw(::V) = SMatrix{2,2}(0.5 * [1.0+im 1.0-im; 1.0-im 1.0+im])
matrix_rep_raw(::Vi) = SMatrix{2,2}(0.5 * [1.0-im 1.0+im; 1.0+im 1.0-im])
matrix_rep_raw(::S)  = SMatrix{2,2}([1.0 0.0; 0.0 im])
matrix_rep_raw(::Si) = SMatrix{2,2}([1.0 0.0; 0.0 -im])
matrix_rep_raw(::T)  = SMatrix{2,2}([1.0 0.0; 0.0 exp(im * π / 4.0)])
matrix_rep_raw(::Ti) = SMatrix{2,2}([1.0 0.0; 0.0 exp(-im * π / 4.0)])
matrix_rep_raw(g::CNot) = SMatrix{4,4}(complex([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 1.0 0.0]))
matrix_rep_raw(g::CY) = SMatrix{4,4}([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 -im; 0.0 0.0 im 0.0])
matrix_rep_raw(g::CZ) = SMatrix{4,4}(complex([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 -1.0]))
matrix_rep_raw(g::CCNot) = SMatrix{8,8}(complex([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0]))
matrix_rep_raw(g::CSwap) = SMatrix{8,8}(complex([1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]))

matrix_rep(g::I) = matrix_rep_raw(g) 
for G in (:Y, :H, :Swap, :CNot, :CY, :CZ, :CCNot, :CSwap, :GPi, :ECR)
    @eval begin
        function matrix_rep(g::$G)
            iseven(g.pow_exponent) && return matrix_rep_raw(I(), qubit_count(g))
            isodd(g.pow_exponent)  && return matrix_rep_raw(g)
            # this path is taken if n == 2.1, for example
            return SMatrix{2^qubit_count($G), 2^qubit_count($G)}(matrix_rep_raw(g)^g.pow_exponent)
        end
    end
end

for (G, G2) in ((:X, :V), (:Z, :S))
    @eval begin
        function matrix_rep(g::$G)
            iseven(g.pow_exponent) && return matrix_rep_raw(I(), qubit_count(g))
            isodd(g.pow_exponent) && return matrix_rep_raw(g)
            g.pow_exponent == 0.5 && return matrix_rep_raw($G2())
            # this path is taken if n == 2.1, for example
            return SMatrix{2, 2}(matrix_rep_raw(g)^g.pow_exponent)
        end
    end
end

function matrix_rep(g::ISwap)
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return SMatrix{4, 4}([1.0 0.0 0.0 0.0; 0.0 0.0 -im 0.0; 0.0 -im 0.0 0.0; 0.0 0.0 0.0 1.0])
    return SMatrix{4, 4}(matrix_rep_raw(g)^n)
end

function matrix_rep(g::PSwap)
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    isodd(n) && return matrix_rep_raw(PSwap(g.angle[1]*n))
    iseven(n) && return SMatrix{4,4}(diagm([1.0, exp(im*g.angle[1]*n), exp(im*g.angle[1]*n), 1.0]))
    return SMatrix{4, 4}(matrix_rep_raw(g)^n)
end

function matrix_rep(g::PRx)
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return matrix_rep_raw(PRx(-g.angle[1], g.angle[2]))
    return SMatrix{2, 2}(matrix_rep_raw(g)^n)
end

for G in (:Rx, :Ry, :Rz, :PhaseShift, :CPhaseShift, :CPhaseShift00, :CPhaseShift01, :CPhaseShift10, :XX, :YY, :ZZ, :XY, :SingleExcitation, :DoubleExcitation, :MultiRZ)
    @eval function matrix_rep(g::$G)
        n = g.pow_exponent
        iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
        isone(n) && return matrix_rep_raw(g)
        isinteger(n) && return matrix_rep_raw($G(g.angle[1]*n))
        return SMatrix{2^qubit_count(g),2^qubit_count(g)}(matrix_rep_raw($G(g.angle[1])) ^ n)
    end
end
    
function matrix_rep(g::CV)
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return SMatrix{4, 4}(complex(ComplexF64[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.5-0.5im 0.5+0.5im; 0.0 0.0 0.5+0.5im 0.5-0.5im]))
    return SMatrix{4, 4}(matrix_rep_raw(g)^n)
end

function matrix_rep(g::MS)
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return matrix_rep_raw(MS(g.angle[1] + π, g.angle[2], g.angle[3]))
    return SMatrix{4, 4}(matrix_rep_raw(g)^n)
end

for (G, Gi, G2) in ((:V, :Vi, :X), (:S, :Si, :Z))
    @eval begin
        function matrix_rep(g::$G)
            n = g.pow_exponent
            mod(n, 4) == 0 && return matrix_rep_raw(I(), qubit_count(g))
            mod(n, 4) == 1 && return matrix_rep_raw($G())
            mod(n, 4) == 2 && return matrix_rep_raw($G2())
            mod(n, 4) == 3 && return matrix_rep_raw($Gi())
            # this path is taken if n == 2.1, for example
            return SMatrix{2, 2}(matrix_rep_raw(g)^n)
        end
        function matrix_rep(g::$Gi)
            n = g.pow_exponent
            mod(n, 4) == 0 && return matrix_rep_raw(I(), qubit_count(g))
            mod(n, 4) == 1 && return matrix_rep_raw($Gi())
            mod(n, 4) == 2 && return matrix_rep_raw($G2())
            mod(n, 4) == 3 && return matrix_rep_raw($G())
            # this path is taken if n == 2.1, for example
            return SMatrix{2, 2}(matrix_rep_raw(g)^n)
        end
    end
end

function matrix_rep_raw(g::U)::SMatrix{2, 2, ComplexF64}
    θ, ϕ, λ = g.angle
    return SMatrix{2, 2, ComplexF64}([cos(θ/2) -exp(im*λ)*sin(θ/2); exp(im*ϕ)*sin(θ/2) exp(im*(λ+ϕ))*cos(θ/2)])
end

function matrix_rep(g::U)::SMatrix{2, 2, ComplexF64}
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return matrix_rep_raw(U(-g.angle[1], -g.angle[3], -g.angle[2]))
    return SMatrix{2, 2}(matrix_rep_raw(g) ^ n)
end
function matrix_rep(g::T)::SMatrix{2, 2, ComplexF64}
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return matrix_rep_raw(Ti())
    return matrix_rep(PhaseShift(n*π/4))
end
function matrix_rep(g::Ti)::SMatrix{2, 2, ComplexF64}
    n = g.pow_exponent
    iszero(n) && return matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == -1 && return matrix_rep_raw(T())
    return matrix_rep(PhaseShift(-n*π/4))
end

matrix_rep_raw(g::CV) = SMatrix{4,4}([1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.5+0.5im 0.5-0.5im; 0.0 0.0 0.5-0.5im 0.5+0.5im])

matrix_rep_raw(g::PRx) = SMatrix{2,2}(
    [
         cos(g.angle[1] / 2.0) -im*exp(-im*g.angle[2])*sin(g.angle[1]/2.0)
        -im*exp(im*g.angle[2])*sin(g.angle[1]/2.0) cos(g.angle[1] / 2.0)
    ],
)
matrix_rep_raw(g::Rz) =
    SMatrix{2,2}([exp(-im * g.angle[1] / 2.0) 0.0; 0.0 exp(im * g.angle[1] / 2.0)])
matrix_rep_raw(g::Rx) = SMatrix{2,2}(
    [
        cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0)
        -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0)
    ],
)
matrix_rep_raw(g::Ry) = SMatrix{2,2}(
    complex(
        [
            cos(g.angle[1] / 2.0) -sin(g.angle[1] / 2.0)
            sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0)
        ],
    ),
)
matrix_rep_raw(g::GPi) =
    SMatrix{2,2}(complex([0 exp(-im * g.angle[1]); exp(im * g.angle[1]) 0]))

matrix_rep_raw(g::GPi2) =
    SMatrix{2,2}(1/√2 * complex([1.0 -im*exp(-im * g.angle[1]); -im*exp(im * g.angle[1]) 1.0]))
matrix_rep_raw(g::MS) = SMatrix{4,4}(
    complex(
        [
            cos(g.angle[3] / 2.0) 0.0 0.0 -im*exp(-im * (g.angle[1] + g.angle[2]))*sin(g.angle[3] / 2);
            0.0  cos(g.angle[3] / 2.0) -im*exp(-im * (g.angle[1] - g.angle[2]))*sin(g.angle[3] / 2.0) 0.0;
            0.0 -im*exp(im * (g.angle[1] - g.angle[2]))*sin(g.angle[3] / 2) cos(g.angle[3] / 2.0) 0.0;
            -im*exp(im * (g.angle[1] + g.angle[2]))*sin(g.angle[3] / 2.0) 0.0 0.0 cos(g.angle[3] / 2.0);
        ],
    ),
)
matrix_rep_raw(g::PhaseShift) = SMatrix{2,2}([1.0 0.0; 0.0 exp(im * g.angle[1])])

# controlled unitaries
for (cph, ind) in
    ((:CPhaseShift, 4), (:CPhaseShift00, 1), (:CPhaseShift01, 3), (:CPhaseShift10, 2))
    @eval begin
        matrix_rep_raw(g::$cph) = Diagonal(
            SVector{4,ComplexF64}(
                setindex!(ones(ComplexF64, 4), exp(im * g.angle[1]), $ind),
            ),
        )
    end
end

for (sw, factor) in ((:Swap, 1.0), (:ISwap, im), (:PSwap, :(exp(im * g.angle[1]))))
    @eval begin
        matrix_rep_raw(g::$sw) = SMatrix{4,4,ComplexF64}(
            [1.0 0.0 0.0 0.0 0.0 0.0 $factor 0.0 0.0 $factor 0.0 0.0 0.0 0.0 0.0 1.0],
        )
    end
end
matrix_rep_raw(g::XX) = SMatrix{4,4}(
    [
        cos(g.angle[1] / 2.0) 0.0 0.0 -im*sin(g.angle[1] / 2)
        0.0 cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0) 0.0
        0.0 -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        -im*sin(g.angle[1] / 2.0) 0.0 0.0 cos(g.angle[1] / 2.0)
    ],
)
matrix_rep_raw(g::XY) = SMatrix{4,4}(
    [
        1.0 0.0 0.0 0.0
        0.0 cos(g.angle[1] / 2.0) im*sin(g.angle[1] / 2.0) 0.0
        0.0 im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        0.0 0.0 0.0 1.0
    ],
)
matrix_rep_raw(g::YY) = SMatrix{4,4}(
    [
        cos(g.angle[1] / 2.0) 0.0 0.0 im*sin(g.angle[1] / 2)
        0.0 cos(g.angle[1] / 2.0) -im*sin(g.angle[1] / 2.0) 0.0
        0.0 -im*sin(g.angle[1] / 2.0) cos(g.angle[1] / 2.0) 0.0
        im*sin(g.angle[1] / 2.0) 0.0 0.0 cos(g.angle[1] / 2.0)
    ],
)
matrix_rep_raw(g::ZZ) = SMatrix{4,4}(
    diagm([
        exp(-im * g.angle[1] / 2.0),
        exp(im * g.angle[1] / 2.0),
        exp(im * g.angle[1] / 2.0),
        exp(-im * g.angle[1] / 2.0),
    ]),
)
# 1/√2 * (IX - XY)
matrix_rep_raw(g::ECR) = SMatrix{4,4}(1/√2 * [0 1 0 im; 1 0 -im 0; 0 im 0 1; -im 0 1 0])
matrix_rep_raw(g::Unitary) = g.matrix
function matrix_rep(g::Gate)
    n = g.pow_exponent
    iszero(n) && matrix_rep_raw(I(), qubit_count(g))
    isone(n) && return matrix_rep_raw(g)
    n == 0.5 && return √matrix_rep_raw(g)
    n == -0.5 && return inv(√matrix_rep_raw(g))
    return SMatrix{2^qubit_count(g),2^qubit_count(g)}( matrix_rep_raw(g)^n )
end

apply_gate!(::Val{false}, g::I, state_vec::AbstractStateVector{T}, qubits::Int...) where {T<:Complex} =
    return
apply_gate!(::Val{true}, g::I, state_vec::AbstractStateVector{T}, qubits::Int...) where {T<:Complex} =
    return
apply_gate!(::Measure, state_vec, args...) = return

function apply_gate!(
    g_matrix::Union{SMatrix{2,2,T}, Diagonal{T,SVector{2,T}}},
    state_vec::AbstractStateVector{T},
    qubit::Int,
) where {T<:Complex}
    n_amps, endian_qubit = get_amps_and_qubits(state_vec, qubit)
    n_tasks  = n_amps >> 1
    # We split the set of amplitude pairs (`n_tasks`) into `n_chunks` groups
    # so that we do not spawn too many tasks at each gate invocation. Each
    # iteration of the threaded for loop spawns a `Task`, which has a (small)
    # startup and teardown time. If the lifetime of the `Task` is short, these
    # overheads can come to dominate the overall runtime.
    # Using this chunking protects the Julia scheduler from being overwhelmed
    # and prevents the time to spawn and reap Julia `Task`s from dominating
    # the runtime.
    n_chunks = max(div(n_tasks, CHUNK_SIZE), 1)
    # Pairs amplitudes for the gate application. For example,
    # 000000 and 000100 would be paired by `flipper` if `endian_qubit` is 3.
    flipper  = 1 << endian_qubit
    # check if the qubit index to flip to find amplitude pairs is larger or smaller
    # than the chunk size.
    is_small_target = flipper < CHUNK_SIZE
    g_00, g_10, g_01, g_11 = g_matrix
    Threads.@threads for chunk_index = 0:n_chunks-1
        # my_amps is the group of amplitude generators which this `Task` will process
        my_amps = if n_chunks > 1
            chunk_index*CHUNK_SIZE:((chunk_index+1)*CHUNK_SIZE-1)
        else
            0:n_tasks-1
        end
        lower_ix  = pad_bit(my_amps[1], endian_qubit) + 1
        higher_ix = lower_ix + flipper
        for task_amp = 0:length(my_amps)-1
            if is_small_target && div(task_amp, flipper) > 0 && mod(task_amp, flipper) == 0
                lower_ix  = higher_ix
                higher_ix = lower_ix + flipper
            end
            @inbounds begin
                lower_amp  = state_vec[lower_ix]
                higher_amp = state_vec[higher_ix]
                state_vec[lower_ix]  = g_00 * lower_amp + g_01 * higher_amp
                state_vec[higher_ix] = g_10 * lower_amp + g_11 * higher_amp
            end
            lower_ix  += 1
            higher_ix += 1
        end
    end
end

# generic two-qubit non controlled unitaries
function apply_gate!(
    g_mat::M,
    state_vec::AbstractStateVector{T},
    t1::Int,
    t2::Int,
) where {T<:Complex,M<:Union{SMatrix{4,4,T},Diagonal{T,SVector{4,T}}}}
    n_amps, (endian_t1, endian_t2) = get_amps_and_qubits(state_vec, t1, t2)
    small_t, big_t = minmax(endian_t1, endian_t2)
    chunk_size     = CHUNK_SIZE
    n_tasks        = n_amps >> 2
    n_chunks       = max(div(n_tasks, CHUNK_SIZE), 1)
    chunked_amps   = collect(Iterators.partition(0:n_tasks-1, CHUNK_SIZE))
    Threads.@threads for chunk_index = 1:n_chunks
        for ix in chunked_amps[chunk_index]
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

function apply_controlled_gate!(
    g_matrix::Union{SMatrix{2,2,T}, Diagonal{T,SVector{2,T}}, Matrix{T}},
    c_bit::Bool,
    state_vec::AbstractStateVector{T},
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
    g_matrix::Union{SMatrix{4, 4, T}, Diagonal{T, SVector{4, T}}, Matrix{T}},
    c_bit::Bool,
    state_vec::AbstractStateVector{T},
    control::Int,
    target_1::Int,
    target_2::Int,
) where {T<:Complex}
    n_amps, (endian_control, endian_t1, endian_t2) =
        get_amps_and_qubits(state_vec, control, target_1, target_2)
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
    g_matrix::Union{SMatrix{2, 2, T}, Diagonal{T, SVector{2, T}}, Matrix{T}},
    c1_bit::Bool,
    c2_bit::Bool,
    state_vec::AbstractStateVector{T},
    control_1::Int,
    control_2::Int,
    target::Int,
) where {T<:Complex}
    n_amps, (endian_c1, endian_c2, endian_target) =
        get_amps_and_qubits(state_vec, control_1, control_2, target)
    small_t, mid_t, big_t = sort([endian_c1, endian_c2, endian_target])
    g_00, g_10, g_01, g_11 = g_matrix
    Threads.@threads for ix = 0:div(n_amps, 8)-1
        # insert 0 at c1, 0 at c2, 0 at target
        padded_ix = pad_bits(ix, [small_t, mid_t, big_t])
        # flip c1 and c2 if needed
        lower_ix = padded_ix
        for (c_bit, c_qubit) in zip((c1_bit, c2_bit), (endian_c1, endian_c2))
            c_bit && (lower_ix = flip_bit(lower_ix, c_qubit))
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
        apply_gate!(::Val{$V}, gate::Control{G, B}, state_vec::AbstractStateVector{T}, qubits::Int...) where {T<:Complex, G<:Gate, B} = apply_controlled_gate!(Val($V), Val(B), gate, gate.g ^ gate.pow_exponent, state_vec, gate.bitvals, qubits...)
        apply_gate!(::Val{$V}, gate::Control{MultiQubitPhaseShift{N}, B}, state_vec::AbstractStateVector{T}, qubits::Int...) where {T<:Complex, N, B} = apply_gate!($f(im), gate, state_vec, qubits...)
        apply_gate!(::Val{$V}, gate::Unitary, state_vec::AbstractStateVector{T}, targets::Vararg{Int,NQ}) where {T<:Complex,NQ} = apply_gate!($f(SMatrix{2^NQ, 2^NQ, ComplexF64}(matrix_rep(gate))), state_vec, targets...)
        apply_gate!(::Val{$V}, g::G, state_vec::AbstractStateVector{T}, qubits::Int...) where {G<:Gate,T<:Complex} = apply_gate!($f(matrix_rep(g)), state_vec, qubits...)
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{1},
            # controlled gate on *all* qubits (e.g. CV)
            gate::G,
            # gate on *target* qubit bits only (e.g. V)
            target_gate::TG,
            state_vec::AbstractStateVector{T},
            control_bits::NTuple{1,Int},
            control::Int,
            target::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(target_gate)), control_bits[1] == 1, state_vec, control, target)
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{1},
            gate::G,
            target_gate::TG,
            state_vec::AbstractStateVector{T},
            control_bits::NTuple{1,Int},
            control::Int,
            target_1::Int,
            target_2::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(target_gate)), control_bits[1] == 1, state_vec, control, target_1, target_2)
        # doubly controlled unitaries
        apply_controlled_gate!(
            ::Val{$V},
            ::Val{2},
            gate::G,
            target_gate::TG,
            state_vec::AbstractStateVector{T},
            control_bits::NTuple{2,Int},
            control_1::Int,
            control_2::Int,
            target::Int,
        ) where {G<:Gate,TG<:Gate,T<:Complex} = apply_controlled_gate!($f(matrix_rep(target_gate)), control_bits[1] == 1, control_bits[2] == 1, state_vec, control_1, control_2, target)
    end
end
apply_gate!(g::G, args...) where {G<:Gate} = apply_gate!(Val(false), g, args...)

for (control_gate, target_gate, n_controls) in (
        (:CNot, :X, 1),
        (:CY, :Y, 1),
        (:CZ, :Z, 1),
        (:CV, :V, 1),
        (:CSwap, :Swap, 1),
        (:CCNot, :X, 2),
    ), Vc in (false, true)

    @eval begin
        apply_gate!(::Val{$Vc}, gate::$control_gate, state_vec::AbstractStateVector{T}, qubits::Int...,) where {T<:Complex} = apply_controlled_gate!(Val($Vc), Val($n_controls), gate, $target_gate() ^ gate.pow_exponent, state_vec, ntuple(control_bit->1, Val($n_controls)), qubits...)
    end
end

function apply_gate!(
    factor::Complex,
    gate::Control{MultiQubitPhaseShift{N}, B},
    state_vec::AbstractStateVector{T},
    qubits::Int...,
) where {T<:Complex, N, B}
    bit_values  = gate.bitvals
    g_matrix    = ones(ComplexF64, 2^N)
    phase_index = 1
    # find index of state to shift phase of
    for (bit_index, bit_value) in enumerate(bit_values)
        phase_index += bit_value << (bit_index - 1)
    end
    g_matrix[phase_index] = exp(factor*gate.g.angle[1]*gate.g.pow_exponent * gate.pow_exponent)
    apply_gate!(Diagonal(SVector{2^N, ComplexF64}(g_matrix)), state_vec, qubits...)
end

function apply_gate!(
    g_matrix::Union{SMatrix{N, N, T}, Diagonal{T, SVector{N, T}}},
    state_vec::AbstractStateVector{T},
    targets::Vararg{Int,NQ},
) where {T<:Complex,NQ,N}
    n_amps, endian_ts = get_amps_and_qubits(state_vec, targets...)
    ordered_ts = sort(collect(endian_ts))
    flip_list  = map(0:N-1) do target_amp
        # select qubits to flip for each index in the target space
        to_flip = Bool[(((1 << qubit_ix) & target_amp) >> qubit_ix) for qubit_ix = 0:NQ-1]
        return ordered_ts[to_flip]
    end
    Threads.@threads for ix = 0:div(n_amps, N)-1
        padded_ix = pad_bits(ix, ordered_ts)
        ixs = SVector{N,Int}(flip_bits(padded_ix, to_flip) + 1 for to_flip in flip_list)
        @views @inbounds begin
            amps = SVector{N,T}(state_vec[ixs])
            new_amps = g_matrix * amps
            state_vec[ixs] = new_amps
        end
    end
    return
end
