function Base.:(^)(g::Gate, n::Real)
    iszero(n) && return I()
    isone(n) && return g
    return Unitary(Matrix(Matrix(matrix_rep(g))^n))
end

Base.:(^)(g::I, n::Real) = I()
for G in (:X, :Y, :Z, :H, :Swap, :CNot, :CY, :CZ, :CCNot, :CSwap, :GPi, :ECR)
    @eval begin
        function Base.:(^)(g::$G, n::Real)
            iseven(n) && return I()
            isodd(n) && return g
            # this path is taken if n == 2.1, for example
            return Unitary(Matrix(matrix_rep(g)^n))
        end
    end
end
for G in (
    :Rx,
    :Ry,
    :Rz,
    :PhaseShift,
    :CPhaseShift,
    :CPhaseShift00,
    :CPhaseShift01,
    :CPhaseShift10,
    :XX,
    :YY,
    :ZZ,
    :XY,
)
    @eval function Base.:(^)(g::$G, n::Real)
        iszero(n) && return I()
        isone(n) && return g
        return $G(g.angle[1]*n)
    end
end

function Base.:(^)(g::PSwap, n::Real)
    iszero(n) && return I()
    isone(n) && return g
    isodd(n) && return PSwap(g.angle[1]*n)
    iseven(n) && return Unitary(diagm([1.0, exp(im*g.angle[1]*n), exp(im*g.angle[1]*n), 1.0]))
    return Unitary(Matrix(matrix_rep(g)^n)) 
end

function Base.:(^)(g::MS, n::Real)
    iszero(n) && return I()
    isone(n) && return g
    isinteger(n) && return MS(g.angle[1], g.angle[2], n*g.angle[3])
    return Unitary(Matrix(matrix_rep(g)^n)) 
end

for (G, Gi, G2) in ((:V, :Vi, :X), (:S, :Si, :Z))
    @eval begin
        function Base.:(^)(g::$G, n::Real)
            mod(n, 4) == 0 && return I()
            mod(n, 4) == 1 && return $G()
            mod(n, 4) == 2 && return $G2()
            mod(n, 4) == 3 && return $Gi()
            # this path is taken if n == 2.1, for example
            return Unitary(Matrix(matrix_rep(g)^n))
        end
        function Base.:(^)(g::$Gi, n::Real)
            mod(n, 4) == 0 && return I()
            mod(n, 4) == 1 && return $Gi()
            mod(n, 4) == 2 && return $G2()
            mod(n, 4) == 3 && return $G()
            # this path is taken if n == 2.1, for example
            return Unitary(Matrix(matrix_rep(g)^n))
        end
    end
end

Base.:(^)(g::MultiQubitPhaseShift{N}, n::Real) where {N} = MultiQubitPhaseShift{N}((g.angle[1]*n,))
function Base.:(^)(g::Control{<:Gate, B}, n::Real) where {B}
    iszero(n) && return Braket.I()
    isone(n) && return g
    return Control(g.g ^ n, g.bitvals)
end
Base.:(^)(g::Control{MultiQubitPhaseShift{N}, B}, n::Real) where {N, B} = Control{MultiQubitPhaseShift{N}, B}(g.g ^ n, g.bitvals)
