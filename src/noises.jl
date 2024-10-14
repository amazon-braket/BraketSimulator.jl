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
StructTypes.constructfrom(::Type{Kraus}, nt::Quasar.CircuitInstruction) = Kraus(nt.arguments)

for (typ, label, arg, qc) in ((:BitFlip, "Bit flip", :probability, 1),
                              (:PhaseFlip, "Phase flip", :probability, 1),
                              (:Depolarizing, "Depolarizing", :probability, 1),
                              (:AmplitudeDamping, "Amplitude damping", :gamma, 1),
                              (:PhaseDamping, "Phase damping", :gamma, 1),
                              (:TwoQubitDepolarizing, "Two qubit depolarizing", :probability, 2),
                              (:TwoQubitDephasing, "Two qubit dephasing", :probability, 2),
                             )
    @eval begin
        @doc """
            $($typ) <: Noise

        $label noise operation.
        """
        struct $typ <: Noise
            $arg::Union{Float64, FreeParameter}
        end
        Parametrizable(::$typ) = Parametrized()
        qubit_count(::$typ) = $qc
        StructTypes.constructfrom(::Type{$typ}, nt::Quasar.CircuitInstruction) = $typ(only(nt.arguments))
    end
end

"""
    PauliChannel <: Noise

PauliChannel noise operation.
"""
struct PauliChannel <: Noise
    probX::Union{Float64, FreeParameter}
    probY::Union{Float64, FreeParameter}
    probZ::Union{Float64, FreeParameter}
end
Parametrizable(::PauliChannel) = Parametrized()
qubit_count(::PauliChannel) = 1
StructTypes.constructfrom(::Type{PauliChannel}, nt::Quasar.CircuitInstruction) = PauliChannel(nt.arguments...)

"""
    GeneralizedAmplitudeDamping <: Noise

GeneralizedAmplitudeDamping noise operation.
"""
struct GeneralizedAmplitudeDamping <: Noise
    probability::Union{Float64, FreeParameter}
    gamma::Union{Float64, FreeParameter}
end
Parametrizable(::GeneralizedAmplitudeDamping) = Parametrized()
qubit_count(::GeneralizedAmplitudeDamping) = 1
StructTypes.constructfrom(::Type{GeneralizedAmplitudeDamping}, nt::Quasar.CircuitInstruction) = GeneralizedAmplitudeDamping(nt.arguments...)

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
qubit_count(::MultiQubitPauliChannel{N}) where {N} = N
Parametrizable(::MultiQubitPauliChannel) = Parametrized()
function MultiQubitPauliChannel(probabilities::Dict{String, <:Union{Float64, FreeParameter}})
    N = length(first(keys(probabilities)))
    return MultiQubitPauliChannel{N}(probabilities)
end
Base.:(==)(c1::MultiQubitPauliChannel{N}, c2::MultiQubitPauliChannel{M}) where {N,M} = N == M && c1.probabilities == c2.probabilities

Parametrizable(::Noise) = NonParametrized()
parameters(g::Noise)     = parameters(Parametrizable(g), g)
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
