"""
    Noise <: QuantumOperator

Abstract type representing a quantum noise operation.
"""
abstract type Noise <: QuantumOperator end

"""
    Kraus <: Noise

Kraus noise operation.
"""
struct Kraus <: Noise
    matrices::Vector{Matrix{ComplexF64}}
end
Base.:(==)(k1::Kraus, k2::Kraus) = k1.matrices == k2.matrices
qubit_count(g::Kraus) = convert(Int, log2(size(g.matrices[1], 1)))

"""
    BitFlip <: Noise

BitFlip noise operation.
"""
struct BitFlip <: Noise
    probability::Union{Float64, FreeParameter}
end
Parametrizable(g::BitFlip) = Parametrized()
qubit_count(g::BitFlip) = 1

"""
    PhaseFlip <: Noise

PhaseFlip noise operation.
"""
struct PhaseFlip <: Noise
    probability::Union{Float64, FreeParameter}
end
Parametrizable(g::PhaseFlip) = Parametrized()
qubit_count(g::PhaseFlip) = 1

"""
    PauliChannel <: Noise

PauliChannel noise operation.
"""
struct PauliChannel <: Noise
    probX::Union{Float64, FreeParameter}
    probY::Union{Float64, FreeParameter}
    probZ::Union{Float64, FreeParameter}
end
Parametrizable(g::PauliChannel) = Parametrized()
qubit_count(g::PauliChannel) = 1

"""
    AmplitudeDamping <: Noise

AmplitudeDamping noise operation.
"""
struct AmplitudeDamping <: Noise
    gamma::Union{Float64, FreeParameter}
end
Parametrizable(g::AmplitudeDamping) = Parametrized()
qubit_count(g::AmplitudeDamping) = 1

"""
    PhaseDamping <: Noise

PhaseDamping noise operation.
"""
struct PhaseDamping <: Noise
    gamma::Union{Float64, FreeParameter}
end
Parametrizable(g::PhaseDamping) = Parametrized()
qubit_count(g::PhaseDamping) = 1

"""
    Depolarizing <: Noise

Depolarizing noise operation.
"""
struct Depolarizing <: Noise
    probability::Union{Float64, FreeParameter}
end
Parametrizable(g::Depolarizing) = Parametrized()
qubit_count(g::Depolarizing) = 1

"""
    TwoQubitDephasing <: Noise

TwoQubitDephasing noise operation.
"""
struct TwoQubitDephasing <: Noise
    probability::Union{Float64, FreeParameter}
end
Parametrizable(g::TwoQubitDephasing) = Parametrized()
qubit_count(g::TwoQubitDephasing) = 2

"""
    TwoQubitDepolarizing <: Noise

TwoQubitDepolarizing noise operation.
"""
struct TwoQubitDepolarizing <: Noise
    probability::Union{Float64, FreeParameter}
end
Parametrizable(g::TwoQubitDepolarizing) = Parametrized()
qubit_count(g::TwoQubitDepolarizing) = 2

"""
    GeneralizedAmplitudeDamping <: Noise

GeneralizedAmplitudeDamping noise operation.
"""
struct GeneralizedAmplitudeDamping <: Noise
    probability::Union{Float64, FreeParameter}
    gamma::Union{Float64, FreeParameter}
end
Parametrizable(g::GeneralizedAmplitudeDamping) = Parametrized()
qubit_count(g::GeneralizedAmplitudeDamping) = 1

"""
    MultiQubitPauliChannel{N} <: Noise

Pauli channel noise operation on `N` qubits.
"""
struct MultiQubitPauliChannel{N} <: Noise
    probabilities::Dict{String, Union{Float64, FreeParameter}}
end
"""
    TwoQubitPauliChannel <: Noise

Pauli channel noise operation on two qubits.
"""
TwoQubitPauliChannel = MultiQubitPauliChannel{2} 
qubit_count(g::MultiQubitPauliChannel{N}) where {N} = N
Parametrizable(g::MultiQubitPauliChannel) = Parametrized()
function MultiQubitPauliChannel(probabilities::Dict{String, <:Union{Float64, FreeParameter}})
    N = length(first(keys(probabilities)))
    return MultiQubitPauliChannel{N}(probabilities)
end
Base.:(==)(c1::MultiQubitPauliChannel{N}, c2::MultiQubitPauliChannel{M}) where {N,M} = N == M && c1.probabilities == c2.probabilities

Parametrizable(g::Noise) = NonParametrized()
parameters(g::Noise) = parameters(Parametrizable(g), g)
parameters(::Parametrized, g::N) where {N<:Noise} = filter(x->x isa FreeParameter, [getproperty(g, fn) for fn in fieldnames(N)])
parameters(::NonParametrized, g::Noise) = FreeParameter[]

# nosemgrep
function bind_value!(::Parametrized, g::N, params::Dict{Symbol, <:Real}) where {N<:Noise}
    new_args = OrderedDict(zip(fieldnames(N), (getproperty(g, fn) for fn in fieldnames(N)))) 
    for fp in findall(v->v isa FreeParameter, new_args)
        new_args[fp] = get(params, getproperty(g, fp).name, new_args[fp])
    end
    return N(values(new_args)...)
end
