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

module IR
import ..BraketSimulator: braketSchemaHeader, BraketSchemaBase, AbstractProgram
using StructTypes

export Program, AbstractIR, AbstractProgramResult, IRObservable

abstract type AbstractIR end

const IRObservable = Union{Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, String}

abstract type AbstractProgramResult <: AbstractIR end
StructTypes.StructType(::Type{AbstractProgramResult}) = StructTypes.AbstractType()
StructTypes.subtypes(::Type{AbstractProgramResult}) = (amplitude=Amplitude, expectation=Expectation, probability=Probability, sample=Sample, statevector=StateVector, densitymatrix=DensityMatrix, variance=Variance, adjoint_gradient=AdjointGradient)

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
end

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
