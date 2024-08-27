using StructTypes

abstract type BraketSchemaBase end

@enum ExponentType int float
ExponentTypeDict = Dict(string(inst)=>inst for inst in instances(ExponentType))

struct braketSchemaHeader
    name::String
    version::String
end

abstract type DeviceActionProperties <: BraketSchemaBase end

abstract type AbstractProgram <: BraketSchemaBase end
StructTypes.StructType(::Type{AbstractProgram}) = StructTypes.AbstractType()
function bsh_f(raw_symbol::Symbol)
    raw = String(raw_symbol)
    v   = eval(Meta.parse(raw))
    return v["name"] == "braket.ir.openqasm.program" ? OpenQasmProgram : Program
end
StructTypes.subtypes(::Type{AbstractProgram}) = StructTypes.SubTypeClosure(bsh_f)
StructTypes.subtypekey(::Type{AbstractProgram}) = :braketSchemaHeader

module IR
import ..BraketSimulator: braketSchemaHeader, BraketSchemaBase, AbstractProgram
using StructTypes

export Program, AbstractIR, AbstractProgramResult, IRObservable

abstract type AbstractIR end

const IRObservable = Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}

abstract type AbstractProgramResult <: AbstractIR end
StructTypes.StructType(::Type{AbstractProgramResult}) = StructTypes.AbstractType()
StructTypes.subtypes(::Type{AbstractProgramResult}) = (probability=Probability, statevector=StateVector, adjoint_gradient=AdjointGradient, densitymatrix=DensityMatrix, sample=Sample, expectation=Expectation, variance=Variance, amplitude=Amplitude)

struct Program <: AbstractProgram
    braketSchemaHeader::braketSchemaHeader
    instructions::Vector{<:Any}
    results::Union{Nothing, Vector{AbstractProgramResult}}
    basis_rotation_instructions::Union{Nothing, Vector{<:Any}}
end

struct Sample <: AbstractProgramResult
    observable::Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}
    targets::Union{Nothing, Vector{Int}}
    type::String
end

struct DensityMatrix <: AbstractProgramResult
    targets::Union{Nothing, Vector{Int}}
    type::String
end

struct Expectation <: AbstractProgramResult
    observable::Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}
    targets::Union{Nothing, Vector{Int}}
    type::String
end

struct Amplitude <: AbstractProgramResult
    states::Vector{String}
    type::String
end

struct Probability <: AbstractProgramResult
    targets::Union{Nothing, Vector{Int}}
    type::String
end

struct AdjointGradient <: AbstractProgramResult
    parameters::Union{Nothing, Vector{String}}
    observable::Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}
    targets::Union{Nothing, Vector{Vector{Int}}}
    type::String
end

struct StateVector <: AbstractProgramResult
    type::String
end

struct Variance <: AbstractProgramResult
    observable::Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}
    targets::Union{Nothing, Vector{Int}}
    type::String
end

end # module

using .IR

@enum DeviceActionType openqasm jaqcd
DeviceActionTypeDict = Dict(string(inst)=>inst for inst in instances(DeviceActionType))

struct ResultType
    name::String
    observables::Union{Nothing, Vector{String}}
    minShots::Union{Nothing, Int}
    maxShots::Union{Nothing, Int}
end

struct AdditionalMetadata
    action::AbstractProgram
    dwaveMetadata::Nothing
    ionqMetadata::Nothing
    rigettiMetadata::Nothing
    oqcMetadata::Nothing
    xanaduMetadata::Nothing
    queraMetadata::Nothing
    simulatorMetadata::Nothing
end

struct ResultTypeValue
    type::AbstractProgramResult
    value::Union{Vector, Float64, Dict}
    ResultTypeValue(type, value::Float64)             = new(type, value)
    ResultTypeValue(type, value::Dict{String, <:Any}) = new(type, value)
    ResultTypeValue(type, value::Vector)              = new(type, value)
end
# for custom lowering of ComplexF64
# just build the tuples directly to
# avoid allocing so much
StructTypes.StructType(::Type{ResultTypeValue}) = StructTypes.CustomStruct()
function StructTypes.lower(rtv::ResultTypeValue)
    lower_complex(v::Complex)            = (real(v), imag(v))
    lower_complex(v::Float64)            = v
    lower_complex(v::Vector{Float64})    = v
    lower_complex(v::Vector{ComplexF64}) = reinterpret(Tuple{Float64, Float64}, v)
    lower_complex(v::Vector)             = map(lower_complex, v)
    lower_complex(v::Dict)               = Dict(k=>lower_complex(v_) for (k,v_) in v)
    return (type=rtv.type, value=lower_complex(rtv.value))
end
ResultTypeValue(nt::@NamedTuple{type::AbstractProgramResult, value::Union{Float64, Dict, Vector}}) = ResultTypeValue(nt.type, nt.value)
StructTypes.lowertype(::Type{ResultTypeValue}) = @NamedTuple{type::AbstractProgramResult, value::Union{Float64, Dict, Vector}} 

struct JaqcdDeviceActionProperties <: DeviceActionProperties
    version::Vector{String}
    actionType::String
    supportedOperations::Union{Nothing, Vector{String}}
    supportedResultTypes::Union{Nothing, Vector{ResultType}}
    disabledQubitRewiringSupported::Union{Nothing, Bool}
end

struct ControlMod
    name::String
    max_qubits::Union{Nothing, Int}
end

struct NegControlMod
    name::String
    max_qubits::Union{Nothing, Int}
end

struct Power
    name::String
    exponent_types::Vector{ExponentType}
end

struct Inverse
    name::String
end

struct OpenQASMDeviceActionProperties <: DeviceActionProperties
    version::Vector{String}
    actionType::String
    supportedOperations::Union{Nothing, Vector{String}}
    supportedModifiers::Vector{Union{ControlMod, NegControlMod, Power, Inverse}}
    supportedPragmas::Vector{String}
    forbiddenPragmas::Vector{String}
    maximumQubitArrays::Union{Nothing, Int}
    maximumClassicalArrays::Union{Nothing, Int}
    forbiddenArrayOperations::Vector{String}
    requiresAllQubitsMeasurement::Bool
    supportPhysicalQubits::Bool
    requiresContiguousQubitIndices::Bool
    supportsPartialVerbatimBox::Bool
    supportsUnassignedMeasurements::Bool
    disabledQubitRewiringSupported::Bool
    supportedResultTypes::Union{Nothing, Vector{ResultType}}
end

struct GateModelParameters <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    qubitCount::Int
    disableQubitRewiring::Bool
end

struct OpenQasmProgram <: AbstractProgram
    braketSchemaHeader::braketSchemaHeader
    source::String
    inputs::Union{Nothing, Dict{String, Union{String, Float64, Int, Vector{Union{String, Float64, Int}}}}}
end

@enum ExecutionDay everyday weekdays weekends monday tuesday wednesday thursday friday saturday sunday
ExecutionDayDict = Dict(string(inst)=>inst for inst in instances(ExecutionDay))

struct DeviceExecutionWindow
    executionDay::Union{ExecutionDay, String}
    windowStartHour::Dates.Time
    windowEndHour::Dates.Time
end

struct DeviceServiceProperties <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    executionWindows::Vector{DeviceExecutionWindow}
    shotsRange::Tuple{Int, Int}
    deviceCost::Nothing
    deviceDocumentation::Nothing
    deviceLocation::Nothing
    updatedAt::Nothing
    getTaskPollIntervalMillis::Nothing
end
StructTypes.defaults(::Type{DeviceServiceProperties}) = Dict{Symbol, Any}(:braketSchemaHeader => braketSchemaHeader("braket.device_schema.device_service_properties", "1"))

struct GateModelSimulatorParadigmProperties <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    qubitCount::Int
end
StructTypes.defaults(::Type{GateModelSimulatorParadigmProperties}) = Dict{Symbol, Any}(:braketSchemaHeader => braketSchemaHeader("braket.device_schema.simulators.gate_model_simulator_paradigm_properties", "1"))

struct GateModelSimulatorDeviceCapabilities <: BraketSchemaBase
    service::DeviceServiceProperties
    action::Dict{Union{DeviceActionType, String}, Union{OpenQASMDeviceActionProperties, JaqcdDeviceActionProperties}}
    deviceParameters::Dict
    braketSchemaHeader::braketSchemaHeader
    paradigm::GateModelSimulatorParadigmProperties
end
StructTypes.defaults(::Type{GateModelSimulatorDeviceCapabilities}) = Dict{Symbol, Any}(:braketSchemaHeader => braketSchemaHeader("braket.device_schema.simulators.gate_model_simulator_device_capabilities", "1"))

struct GateModelSimulatorDeviceParameters <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    paradigmParameters::GateModelParameters
end

struct TaskMetadata <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    id::String
    shots::Int
    deviceId::String
    deviceParameters::Nothing
    createdAt::Union{Nothing, String}
    endedAt::Union{Nothing, String}
    status::Union{Nothing, String}
    failureReason::Union{Nothing, String}
end

struct GateModelTaskResult <: BraketSchemaBase
    braketSchemaHeader::braketSchemaHeader
    measurements::Union{Nothing, Vector{Vector{Int}}}
    measurementProbabilities::Union{Nothing, Dict{String, Float64}}
    resultTypes::Union{Nothing, Vector{ResultTypeValue}}
    measuredQubits::Union{Nothing, Vector{Int}}
    taskMetadata::TaskMetadata
    additionalMetadata::AdditionalMetadata
end

struct GenericDeviceActionProperties <: DeviceActionProperties
    version::Vector{String}
    actionType::Union{DeviceActionType, String}
end
