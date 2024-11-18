module BraketSimulatorBraketExt

using PrecompileTools

@recompile_invalidations begin
    using BraketSimulator, Braket
end

Braket.name(d::BraketSimulator.AbstractSimulator) = BraketSimulator.name(d)
Braket.properties(d::BraketSimulator.AbstractSimulator) = BraketSimulator.properties(d)
Braket.simulate(d::BraketSimulator.AbstractSimulator, program::Braket.OpenQasmProgram, args...; kwargs...) = convert(Braket.GateModelTaskResult, BraketSimulator.simulate(d, convert(BraketSimulator.OpenQasmProgram, program), args...; kwargs...))
Braket.simulate(d::BraketSimulator.AbstractSimulator, program::Braket.Program, qubit_count::Int, shots::Int; kwargs...) = convert(Braket.GateModelTaskResult, BraketSimulator.simulate(d, convert(BraketSimulator.Program, program), shots; kwargs...))

Base.convert(::Type{Braket.TaskMetadata}, tm::BraketSimulator.TaskMetadata) = Braket.TaskMetadata(Braket.braketSchemaHeader("braket.task_result.task_metadata", "1"), tm.id, tm.shots, tm.deviceId, tm.deviceParameters, tm.createdAt, tm.endedAt, tm.status, tm.failureReason)
Base.convert(::Type{Braket.AdditionalMetadata}, am::BraketSimulator.AdditionalMetadata) = Braket.AdditionalMetadata(convert(Braket.AbstractProgram, am.action), nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
Base.convert(::Type{Braket.ResultTypeValue}, rtv::BraketSimulator.ResultTypeValue) = Braket.ResultTypeValue(convert(Braket.AbstractProgramResult, rtv.type), rtv.value)

# nosemgrep
function Base.convert(::Type{Braket.GateModelTaskResult}, r::BraketSimulator.GateModelTaskResult)
    task_meta = convert(Braket.TaskMetadata, r.taskMetadata)
    addl_meta = convert(Braket.AdditionalMetadata, r.additionalMetadata)
    rts = !isnothing(r.resultTypes) ? Braket.ResultTypeValue[convert(Braket.ResultTypeValue, rt) for rt in r.resultTypes] : nothing

    return Braket.GateModelTaskResult(Braket.braketSchemaHeader("braket.task_result.gate_model_task_result", "1"),
                                      r.measurements,
                                      r.measurementProbabilities,
                                      rts,
                                      r.measuredQubits,
                                      task_meta,
                                      addl_meta)
end

Braket.qubit_count(o::O) where {O<:BraketSimulator.Operator}               = BraketSimulator.qubit_count(o)
Braket.qubit_count(o::O) where {O<:BraketSimulator.Observables.Observable} = BraketSimulator.qubit_count(o)
Braket.qubit_count(::Type{O}) where {O<:BraketSimulator.Operator}          = BraketSimulator.qubit_count(O)

for (braket_sim, simulator_sym) in ((:(Braket.X), :(BraketSimulator.X)),
                                    (:(Braket.Y), :(BraketSimulator.Y)),
                                    (:(Braket.Z), :(BraketSimulator.Z)),
                                    (:(Braket.I), :(BraketSimulator.I)),
                                    (:(Braket.H), :(BraketSimulator.H)),
                                    (:(Braket.V), :(BraketSimulator.V)),
                                    (:(Braket.Vi), :(BraketSimulator.Vi)),
                                    (:(Braket.S), :(BraketSimulator.S)),
                                    (:(Braket.Si), :(BraketSimulator.Si)),
                                    (:(Braket.T), :(BraketSimulator.T)),
                                    (:(Braket.Ti), :(BraketSimulator.Ti)),
                                    (:(Braket.ECR), :(BraketSimulator.ECR)),
                                    (:(Braket.CNot), :(BraketSimulator.CNot)),
                                    (:(Braket.CY), :(BraketSimulator.CY)),
                                    (:(Braket.CZ), :(BraketSimulator.CZ)),
                                    (:(Braket.CV), :(BraketSimulator.CV)),
                                    (:(Braket.CCNot), :(BraketSimulator.CCNot)),
                                    (:(Braket.Swap), :(BraketSimulator.Swap)),
                                    (:(Braket.ISwap), :(BraketSimulator.ISwap)),
                                    (:(Braket.CSwap), :(BraketSimulator.CSwap)),
                                   )
    @eval Base.convert(::Type{BraketSimulator.Operator}, g::$braket_sim) = $simulator_sym()
    @eval Base.convert(::Type{Braket.Operator}, g::$simulator_sym) = $braket_sim()
end

for (braket_sim, simulator_sym) in ((:(Braket.Rx), :(BraketSimulator.Rx)),
                                    (:(Braket.Ry), :(BraketSimulator.Ry)),
                                    (:(Braket.Rz), :(BraketSimulator.Rz)),
                                    (:(Braket.GPi), :(BraketSimulator.GPi)),
                                    (:(Braket.GPi2), :(BraketSimulator.GPi2)),
                                    (:(Braket.XX), :(BraketSimulator.XX)),
                                    (:(Braket.XY), :(BraketSimulator.XY)),
                                    (:(Braket.YY), :(BraketSimulator.YY)),
                                    (:(Braket.ZZ), :(BraketSimulator.ZZ)),
                                    (:(Braket.PSwap), :(BraketSimulator.PSwap)),
                                    (:(Braket.PhaseShift), :(BraketSimulator.PhaseShift)),
                                    (:(Braket.CPhaseShift), :(BraketSimulator.CPhaseShift)),
                                    (:(Braket.CPhaseShift00), :(BraketSimulator.CPhaseShift00)),
                                    (:(Braket.CPhaseShift01), :(BraketSimulator.CPhaseShift01)),
                                    (:(Braket.CPhaseShift10), :(BraketSimulator.CPhaseShift10)),
                                   )
    @eval Base.convert(::Type{BraketSimulator.Operator}, g::$braket_sim) = $simulator_sym(convert(Union{BraketSimulator.FreeParameter, Real}, g.angle[1]))
    @eval Base.convert(::Type{Braket.Operator}, g::$simulator_sym) = $braket_sim(convert(Union{Braket.FreeParameter, Real}, g.angle[1]))
end
for (braket_sim, simulator_sym) in ((:(Braket.BitFlip), :(BraketSimulator.BitFlip)),
                                    (:(Braket.PhaseFlip), :(BraketSimulator.PhaseFlip)),
                                    (:(Braket.Depolarizing), :(BraketSimulator.Depolarizing)),
                                    (:(Braket.TwoQubitDephasing), :(BraketSimulator.TwoQubitDephasing)),
                                    (:(Braket.TwoQubitDepolarizing), :(BraketSimulator.TwoQubitDepolarizing)),
                                   )
    @eval Base.convert(::Type{BraketSimulator.Operator}, g::$braket_sim) = $simulator_sym(convert(Union{BraketSimulator.FreeParameter, Real}, g.probability))
    @eval Base.convert(::Type{Braket.Operator}, g::$simulator_sym) = $braket_sim(convert(Union{Braket.FreeParameter, Real}, g.probability))
end
for (braket_sim, simulator_sym) in ((:(Braket.IR.Expectation), :(BraketSimulator.IR.Expectation)),
                                    (:(Braket.IR.Variance), :(BraketSimulator.IR.Variance)),
                                    (:(Braket.IR.Sample), :(BraketSimulator.IR.Sample)),
                                   )
    @eval Base.convert(::Type{Braket.AbstractProgramResult}, rt::$simulator_sym) = $braket_sim(rt.observable, rt.targets, rt.type)
    @eval Base.convert(::Type{BraketSimulator.AbstractProgramResult}, rt::$braket_sim) = $simulator_sym(rt.observable, rt.targets, rt.type)
end
for (braket_sim, simulator_sym) in ((:(Braket.IR.DensityMatrix), :(BraketSimulator.IR.DensityMatrix)),
                                    (:(Braket.IR.Probability), :(BraketSimulator.IR.Probability)),
                                   )
    @eval Base.convert(::Type{Braket.AbstractProgramResult}, rt::$simulator_sym) = $braket_sim(rt.targets, rt.type)
    @eval Base.convert(::Type{BraketSimulator.AbstractProgramResult}, rt::$braket_sim) = $simulator_sym(rt.targets, rt.type)
end
Base.convert(::Type{Braket.AbstractProgramResult}, rt::BraketSimulator.IR.StateVector) = Braket.IR.StateVector(rt.type)
Base.convert(::Type{BraketSimulator.AbstractProgramResult}, rt::Braket.IR.StateVector) = BraketSimulator.IR.StateVector(rt.type)
Base.convert(::Type{Braket.AbstractProgramResult}, rt::BraketSimulator.IR.Amplitude) = Braket.IR.Amplitude(rt.states, rt.type)
Base.convert(::Type{BraketSimulator.AbstractProgramResult}, rt::Braket.IR.Amplitude) = BraketSimulator.IR.Amplitude(rt.states, rt.type)
Base.convert(::Type{BraketSimulator.Operator}, m::Braket.Measure) = BraketSimulator.Measure(m.index)
Base.convert(::Type{Braket.Operator}, m::BraketSimulator.Measure) = Braket.Measure(m.index)
Base.convert(::Type{BraketSimulator.Operator}, ::Braket.Reset) = BraketSimulator.Reset()
Base.convert(::Type{Braket.Operator}, ::BraketSimulator.Reset) = Braket.Reset()
Base.convert(::Type{BraketSimulator.Operator}, ::Braket.Barrier) = BraketSimulator.Barrier()
Base.convert(::Type{Braket.Operator}, ::BraketSimulator.Barrier) = Braket.Barrier()
Base.convert(::Type{BraketSimulator.Operator}, d::Braket.Delay) = BraketSimulator.Delay(d.duration)
Base.convert(::Type{Braket.Operator}, d::BraketSimulator.Delay) = Braket.Delay(d.duration)

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.MS) = BraketSimulator.MS(ntuple(i->convert(Union{FreeParameter, Real}, g.angle[i]), 3))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.MS) = Braket.MS(ntuple(i->convert(Union{Braket.FreeParameter, Real}, g.angle[i]), 3))

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.U) = BraketSimulator.U(ntuple(i->convert(Union{FreeParameter, Real}, g.angle[i]), 3))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.U) = Braket.U(ntuple(i->convert(Union{Braket.FreeParameter, Real}, g.angle[i]), 3))

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.PRx) = BraketSimulator.PRx(ntuple(i->convert(Union{FreeParameter, Real}, g.angle[i]), 2))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.PRx) = Braket.PRx(ntuple(i->convert(Union{Braket.FreeParameter, Real}, g.angle[i]), 2))

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.Unitary) = BraketSimulator.Unitary(g.matrix)
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.Unitary) = Braket.Unitary(g.matrix)

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.Kraus) = BraketSimulator.Kraus(g.matrices)
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.Kraus) = Braket.Kraus(g.matrices)

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.PauliChannel) = BraketSimulator.PauliChannel(convert(Union{BraketSimulator.FreeParameter,Real}, g.probX), convert(Union{BraketSimulator.FreeParameter,Real}, g.probY), convert(Union{BraketSimulator.FreeParameter,Real}, g.probZ))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.PauliChannel) = Braket.PauliChannel(convert(Union{Braket.FreeParameter,Real}, g.probX), convert(Union{Braket.FreeParameter,Real}, g.probY), convert(Union{Braket.FreeParameter,Real}, g.probZ))

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.GeneralizedAmplitudeDamping) = BraketSimulator.GeneralizedAmplitudeDamping(convert(Union{BraketSimulator.FreeParameter,Real}, g.probability), convert(Union{BraketSimulator.FreeParameter,Real}, g.gamma))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.GeneralizedAmplitudeDamping) = Braket.GeneralizedAmplitudeDamping(convert(Union{Braket.FreeParameter,Real}, g.probability), convert(Union{Braket.FreeParameter,Real}, g.gamma))

Base.convert(::Type{BraketSimulator.Operator}, g::Braket.MultiQubitPauliChannel{N}) where {N} = BraketSimulator.MultiQubitPauliChannel{N}(Dict(k=>convert(Union{BraketSimulator.FreeParameter,Float64}, v) for (k,v) in g.probabilities))
Base.convert(::Type{Braket.Operator}, g::BraketSimulator.MultiQubitPauliChannel{N}) where {N} = Braket.MultiQubitPauliChannel{N}(Dict(k=>convert(Union{Braket.FreeParameter,Real}, v) for (k,v) in g.probabilities))

Base.convert(::Type{BraketSimulator.Instruction}, ix::Braket.Instruction) = BraketSimulator.Instruction(convert(BraketSimulator.Operator, ix.operator), convert(BraketSimulator.QubitSet, ix.target))
Base.convert(::Type{Braket.Instruction}, ix::BraketSimulator.Instruction) = Braket.Instruction(convert(Braket.Operator, ix.operator), convert(Braket.QubitSet, ix.target))

Base.convert(::Type{BraketSimulator.OpenQasmProgram}, p::Braket.OpenQasmProgram) = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), p.source, p.inputs)

# have to handle the special case of 2-qubit Hermitians carefully due to endianness
# nosemgrep
function Base.convert(::Type{BraketSimulator.Program}, p::Braket.Program)
    ixs = [convert(BraketSimulator.Instruction, ix) for ix in p.instructions]
    rts = [convert(BraketSimulator.AbstractProgramResult, rt) for rt in p.results]
    bris = map(p.basis_rotation_instructions) do bri
        if bri.operator isa Unitary && length(bri.target) == 2
            return BraketSimulator.Instruction(BraketSimulator.Unitary(BraketSimulator.fix_endianness(bri.operator.matrix)), bri.target)
        else
            return convert(BraketSimulator.Instruction, bri)
        end
    end
    p = BraketSimulator.Program(BraketSimulator.braketSchemaHeader("braket.ir.jaqcd.program", "1"), ixs, rts, bris)
    return p
end
Base.convert(::Type{Braket.AbstractProgram}, p::BraketSimulator.Program) = Braket.Program(Braket.braketSchemaHeader("braket.ir.jaqcd.program", "1"), [convert(Braket.Instruction, ix) for ix in p.instructions],[convert(Braket.AbstractProgramResult, rt) for rt in p.results], [convert(Braket.Instruction, ix) for ix in p.basis_rotation_instructions])
Base.convert(::Type{Braket.AbstractProgram}, p::BraketSimulator.OpenQasmProgram) = Braket.OpenQasmProgram(Braket.braketSchemaHeader("braket.ir.openqasm.program", "1"), p.source, p.inputs)

function Base.convert(::Type{Braket.Circuit}, c::BraketSimulator.Circuit) 
    ixs = [convert(Braket.Instruction, ix) for ix in c.instructions]
    rts = isempty(c.result_types) ? Braket.Result[] : [convert(Braket.Result, rt) for rt in c.result_types]
    brs = [convert(Braket.Instruction, ix) for ix in c.basis_rotation_instructions] 
    Braket.Circuit(Braket.Moments(), ixs, rts, brs)
end

function __init__()
    Braket._simulator_devices[]["braket_dm_v2"] =
        DensityMatrixSimulator{ComplexF64,Matrix{ComplexF64}}
    Braket._simulator_devices[]["braket_sv_v2"] =
        StateVectorSimulator{ComplexF64,Vector{ComplexF64}}
end

end
