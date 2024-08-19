Base.:(^)(@nospecialize(g::G), n::Real) where {G<:Gate} = G(g.pow_exponent * n)::G
Base.inv(@nospecialize(g::G)) where {G<:Gate} = G(-g.pow_exponent)::G
Base.:(^)(@nospecialize(g::G), n::Real) where {G<:AngledGate} = G(g.angle, g.pow_exponent * n)::G
Base.inv(@nospecialize(g::G)) where {G<:AngledGate} = G(g.angle, -g.pow_exponent)::G
Base.:(^)(g::Unitary, n::Real) = Unitary(g.matrix, g.pow_exponent * n)::Unitary
Base.inv(g::Unitary) = Unitary(g.matrix, -g.pow_exponent)::Unitary
Base.:(^)(@nospecialize(g::Control{<:Gate, B}), n::Real) where {B} = Control(g.g ^ n, g.bitvals)
Base.inv(@nospecialize(g::Control{<:Gate, B})) where {B} = Control(inv(g.g), g.bitvals)