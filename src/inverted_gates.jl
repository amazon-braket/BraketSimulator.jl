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
Base.inv(g::PRx)  = PRx(-g.angle[1], g.angle[2])
Base.inv(g::GPi2) = GPi2(g.angle[1] + π)
Base.inv(g::MS) = MS(g.angle[1] + π, g.angle[2], g.angle[3])

const iswap_inv = ComplexF64[1.0 0.0 0.0 0.0;
                             0.0 0.0 -im 0.0;
                             0.0 -im 0.0 0.0;
                             0.0 0.0 0.0 1.0]
Base.inv(::ISwap) = Unitary(iswap_inv)
const cv_inv = ComplexF64[1.0 0.0 0.0 0.0;
                          0.0 1.0 0.0 0.0;
                          0.0 0.0 0.5-0.5im 0.5+0.5im;
                          0.0 0.0 0.5+0.5im 0.5-0.5im]
Base.inv(::CV) = Unitary(cv_inv)
