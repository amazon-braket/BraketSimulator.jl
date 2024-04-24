for G in (:X, :Y, :Z, :H, :I, :Swap, :CNot, :CY, :CZ, :CCNot, :CSwap, :GPi, :ECR)
    @eval Base.inv(g::$G) = g
end
for G in (
    :Rx,
    :Ry,
    :Rz,
    :PhaseShift,
    :PSwap,
    :CPhaseShift,
    :CPhaseShift00,
    :CPhaseShift01,
    :CPhaseShift10,
    :XX,
    :YY,
    :ZZ,
    :XY,
)
    @eval Base.inv(g::$G) = $G(-g.angle[1])
end
for (G, Gi) in ((:V, :Vi), (:S, :Si), (:T, :Ti))
    @eval begin
        Base.inv(g::$G)  = $Gi()
        Base.inv(g::$Gi) = $G()
    end
end
Base.inv(g::Unitary) = Unitary(inv(g.matrix))
Base.inv(g::GPi2) = GPi2(g.angle[1] + π)
Base.inv(g::MS) = MS(g.angle[1] + π, g.angle[2], g.angle[3])

function Base.inv(g::ISwap)
    u_mat = zeros(ComplexF64, 4, 4)
    u_mat[1, 1] = 1.0
    u_mat[4, 4] = 1.0
    u_mat[2, 3] = -im
    u_mat[3, 2] = -im
    return Unitary(u_mat)
end
function Base.inv(g::CV)
    u_mat = zeros(ComplexF64, 4, 4)
    u_mat[1, 1] = 1.0
    u_mat[2, 2] = 1.0
    u_mat[3, 3] = 0.5-0.5im
    u_mat[3, 4] = 0.5+0.5im
    u_mat[4, 3] = 0.5+0.5im
    u_mat[4, 4] = 0.5-0.5im
    return Unitary(u_mat)
end
