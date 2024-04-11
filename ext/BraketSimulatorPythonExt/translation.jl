convert_ir_matrix(m) = [[[pyconvert(Float64, m__[0]), pyconvert(Float64, m__[1])] for m__ in m_] for m_ in m]

function convert_ir_obs(o)
    if pyisinstance(o, pybuiltins.str)
        return pyconvert(String, o)
    elseif pyisinstance(o, pybuiltins.list)
        jl_o = map(o) do o_
            if pyisinstance(o_, pybuiltins.str)
                return pyconvert(String, o_)
            else
                return convert_ir_matrix(o_)
            end
        end
        return convert(Vector{Union{String, Vector{Vector{Vector{Float64}}}}}, jl_o)
    else
        return convert_ir_matrix(o)
    end
end
function convert_ir_target(t)
    pyis(t, pybuiltins.None) && return nothing
    return [pyconvert(Int, t_) for t_ in t]
end
for (rt, fn) in ((:(Braket.IR.Sample), :jl_convert_sample),
                 (:(Braket.IR.Expectation), :jl_convert_expectation),
                 (:(Braket.IR.Variance), :jl_convert_variance),
                )
    @eval begin
        function $fn(::Type{$rt}, x::Py)
            jl_obs = convert_ir_obs(x.observable)
            jl_ts  = convert_ir_target(x.targets)
            jl_rt  = $rt(jl_obs, jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
        function $fn(::Type{AbstractProgramResult}, x::Py)
            jl_obs = convert_ir_obs(x.observable)
            jl_ts = convert_ir_target(x.targets)
            jl_rt = $rt(jl_obs, jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
    end
end

for (rt, fn) in ((:(Braket.IR.Probability), :jl_convert_probability),
                 (:(Braket.IR.DensityMatrix), :jl_convert_densitymatrix),
                )
    @eval begin
        function $fn(::Type{$rt}, x::Py)
            jl_ts = convert_ir_target(x.targets)
            jl_rt = $rt(jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
        function $fn(::Type{AbstractProgramResult}, x::Py)
            jl_ts = convert_ir_target(x.target)
            jl_rt = $rt(jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
    end
end

function inputs_to_jl(x)
    if pyis(x, pybuiltins.None)
        return nothing
    else
        jl_inputs = map(x.items()) do x_
            k, v = x_
            jl_k = pyconvert(String, k)
            jl_v = pyisinstance(v, pybuiltins.int) ? pyconvert(Int, v) : pyconvert(Float64, v)
            return jl_k=>jl_v
        end
        return Dict(jl_inputs)
    end
end
function jl_convert_oqprogram(::Type{OpenQasmProgram}, x::Py)
    bsh = pyconvert(Braket.braketSchemaHeader, x.braketSchemaHeader)
    source = pyconvert(String, x.source)
    inputs = inputs_to_jl(x.inputs) 
    PythonCall.pyconvert_return(OpenQasmProgram(bsh, source, inputs))
end

function jl_convert_amplitude(::Type{Braket.IR.Amplitude}, x::Py)
    jl_rt = Braket.IR.Amplitude([pyconvert(String, s) for s in x.states], pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_amplitude(::Type{AbstractProgramResult}, x::Py)
    jl_rt = Braket.IR.Amplitude([pyconvert(String, s) for s in x.states], pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end

function jl_convert_statevector(::Type{Braket.IR.StateVector}, x::Py)
    jl_rt = Braket.IR.StateVector(pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_statevector(::Type{AbstractProgramResult}, x::Py)
    jl_rt = Braket.IR.StateVector(pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end

function jl_convert_adjointgradient(::Type{Braket.IR.AdjointGradient}, x::Py)
    jl_targets = pyis(x.targets, pybuiltins.None) ? nothing : [[pyconvert(Int, t__) for t__ in t_] for t_ in t]
    jl_params = pyis(x.targets, pybuiltins.None) ? nothing : [pyconvert(String, p) for p in x.parameters]
    jl_obs = pyisinstance(x.observable, pybuiltins.str) ? pyconvert(String, x.observable) : [convert_ir_obs(o) for o in x.observable] 
    jl_rt = Braket.IR.AdjointGradient(jl_params, jl_obs, jl_targets, pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_adjointgradient(::Type{AbstractProgramResult}, x::Py)
    jl_targets = pyis(x.targets, pybuiltins.None) ? nothing : [[pyconvert(Int, t__) for t__ in t_] for t_ in t]
    jl_params = pyis(x.targets, pybuiltins.None) ? nothing : [pyconvert(String, p) for p in x.parameters]
    jl_obs = pyisinstance(x.observable, pybuiltins.str) ? pyconvert(String, x.observable) : [convert_ir_obs(o) for o in x.observable] 
    jl_rt = Braket.IR.AdjointGradient(jl_params, jl_obs, jl_targets, pyconvert(String, x.type))
end

jl_convert_kraus(::Type{Kraus}, x::Py)::Kraus = PythonCall.pyconvert_return(Kraus([convert_ir_matrix(m) for m in x.matrices]))
jl_convert_kraus(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Kraus([convert_ir_matrix(m) for m in x.matrices]), [pyconvert(Int, t) for t in x.targets]))

jl_convert_unitary(::Type{Unitary}, x::Py)::Unitary = PythonCall.pyconvert_return(Unitary(convert_ir_matrix(x.matrix)))
jl_convert_unitary(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Unitary(convert_ir_matrix(x.matrix)), [pyconvert(Int, t) for t in x.targets]))

jl_convert_cswap(::Type{CSwap}, x::Py)::CSwap = PythonCall.pyconvert_return(CSwap())
jl_convert_cswap(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(CSwap(), vcat(pyconvert(Int, x.control), [pyconvert(Int, t) for t in x.targets])))
jl_convert_ccnot(::Type{CCNot}, x::Py)::CCNot = PythonCall.pyconvert_return(CCNot())
jl_convert_ccnot(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(CCNot(), vcat([pyconvert(Int, c) for c in x.controls], pyconvert(Int, x.target))))

jl_convert_paulichannel(::Type{PauliChannel}, x::Py)::PauliChannel = PythonCall.pyconvert_return(PauliChannel(pyconvert(Float64, x.probX), pyconvert(Float64, x.probY), pyconvert(Float64, x.probZ)))
jl_convert_paulichannel(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(PauliChannel(pyconvert(Float64, x.probX), pyconvert(Float64, x.probY), pyconvert(Float64, x.probZ)), pyconvert(Int, x.target)))

jl_convert_generalizedampdamp(::Type{GeneralizedAmplitudeDamping}, x::Py)::GeneralizedAmplitudeDamping = PythonCall.pyconvert_return(GeneralizedAmplitudeDamping(pyconvert(Float64, x.probability), pyconvert(Float64, x.gamma)))
jl_convert_generalizedampdamp(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(GeneralizedAmplitudeDamping(pyconvert(Float64, x.probability), pyconvert(Float64, x.gamma)), pyconvert(Int, x.target)))

jl_convert_multi_qubit_pauli_channel(::Type{MultiQubitPauliChannel}, x::Py)::MultiQubitPauliChannel = PythonCall.pyconvert_return(MultiQubitPauliChannel(Dict(pyconvert(String, k)=>pyconvert(Float64, v) for (k,v) in x.probabilities)))
jl_convert_multi_qubit_pauli_channel(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(MultiQubitPauliChannel(Dict(pyconvert(String, k)=>pyconvert(Float64, v) for (k,v) in x.probabilities)), x.target))

for (g, fn) in ((:CNot, :jl_convert_cnot),
                (:CY, :jl_convert_cy),
                (:CZ, :jl_convert_cz),
                (:CV, :jl_convert_cv))
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g())
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(), [pyconvert(Int, x.control), pyconvert(Int, x.target)]))
    end
end
for (g, fn) in ((:BitFlip, :jl_convert_bitflip),
                (:PhaseFlip, :jl_convert_phaseflip),
                (:Depolarizing, :jl_convert_depo),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{Float64,FreeParameter}, x.probability)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64,FreeParameter}, x.probability)), pyconvert(Int, x.target)))
    end
end
for (g, fn) in ((:AmplitudeDamping, :jl_convert_amplitudedamp),
                (:PhaseDamping, :jl_convert_phasedamp),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{Float64,FreeParameter}, x.gamma)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64, FreeParameter}, x.gamma)), pyconvert(Int, x.target)))
    end
end
for (g, fn) in ((:TwoQubitDepolarizing, :jl_convert_twoqubitdepo),
                (:TwoQubitDephasing, :jl_convert_twoqubitdeph),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{Float64,FreeParameter}, x.probability)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64, FreeParameter}, x.probability)), [pyconvert(Int, t) for t in x.targets]))
    end
end
for (g, fn) in ((:I, :jl_convert_i),
                (:X, :jl_convert_x),
                (:Y, :jl_convert_y),
                (:Z, :jl_convert_z),
                (:H, :jl_convert_h),
                (:V, :jl_convert_v),
                (:Vi, :jl_convert_vi),
                (:T, :jl_convert_t),
                (:Ti, :jl_convert_ti),
                (:S, :jl_convert_s),
                (:Si, :jl_convert_si)
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g())
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(), pyconvert(Int, x.target)))
    end
end
for (g, fn) in ((:Rx, :jl_convert_rx),
                (:Ry, :jl_convert_ry),
                (:Rz, :jl_convert_rz),
                (:PhaseShift, :jl_convert_phaseshift),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(x.angle))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64, FreeParameter}, x.angle)), pyconvert(Int, x.target)))
    end
end

for (g, fn) in ((:ECR, :jl_convert_ecr),
                (:Swap, :jl_convert_swap),
                (:ISwap, :jl_convert_iswap),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g())
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(), [pyconvert(Int, t) for t in x.targets]))
    end
end

for (g, fn) in ((:XX, :jl_convert_xx),
                (:XY, :jl_convert_xy),
                (:YY, :jl_convert_yy),
                (:ZZ, :jl_convert_zz),
                (:PSwap, :jl_convert_pswap),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(x.angle))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{FreeParameter, Float64}, x.angle)), [pyconvert(Int, t) for t in x.targets]))
    end
end

for (g, fn) in ((:StartVerbatimBox, :jl_convert_startverbatim),
                (:EndVerbatimBox, :jl_convert_endverbatim),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g())
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(), Int[]))
    end
end

for (g, fn) in ((:CPhaseShift, :jl_convert_cphaseshift),
                (:CPhaseShift00, :jl_convert_cphaseshift00),
                (:CPhaseShift10, :jl_convert_cphaseshift10),
                (:CPhaseShift01, :jl_convert_cphaseshift01),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{Float64, FreeParameter}, x.angle)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64, FreeParameter}, x.angle)), [pyconvert(Int, x.control), pyconvert(Int, x.target)]))
    end
end


function jl_convert_ix(::Type{Instruction}, x::Py)::Instruction
    # extract targets
    full_targets = QubitSet()
    for attr ∈ ("control", "target") # single-qubit
        if pyhasattr(x, attr)
            push!(full_targets, pyconvert(Int, pygetattr(x, attr)))
        end
    end   
    for attr ∈ ("controls", "targets") # multi-qubit
        if pyhasattr(x, attr)
            qs = pyconvert(Vector{Int}, pygetattr(x, attr))
            foreach(q->push!(full_targets, q), qs)
        end
    end
    # now build the operator
    sts    = merge(Braket.StructTypes.subtypes(Gate), Braket.StructTypes.subtypes(Noise), Braket.StructTypes.subtypes(CompilerDirective))
    py_name = pyconvert(String, pygetattr(x, "type"))
    jl_type = sts[Symbol(py_name)]
    # bad, can we come up with something more generic
    ix = Instruction(PythonCall.pyconvert(jl_type, x), full_targets)
    PythonCall.pyconvert_return(ix)
end

function jl_convert_bsh(::Type{Braket.braketSchemaHeader}, x::Py)
    name    = pyconvert(String, x.name)
    version = pyconvert(String, x.version)
    PythonCall.pyconvert_return(Braket.braketSchemaHeader(name, version))
end

function jl_convert_circuit(::Type{Braket.IR.Program}, x)
    x_jaqcd      = x._to_jaqcd()
    instructions = [pyconvert(Instruction, ix) for ix in x_jaqcd.instructions]
    results      = pyisinstance(x_jaqcd.results, PythonCall.pybuiltins.list) ? [pyconvert(AbstractProgramResult, rt) for rt in x_jaqcd.results] : AbstractProgramResult[]
    bris         = pyisinstance(x_jaqcd.basis_rotation_instructions, PythonCall.pybuiltins.list) ? [pyconvert(Instruction, ix) for ix in x_jaqcd.basis_rotation_instructions] : Instruction[]
    prog         = Braket.Program(Braket.header_dict[Braket.Program], instructions, results, bris)
    PythonCall.pyconvert_return(prog)
end
function jl_convert_program(::Type{Braket.IR.Program}, x_jaqcd)
    instructions = [pyconvert(Instruction, ix) for ix in x_jaqcd.instructions]
    results      = pyisinstance(x_jaqcd.results, PythonCall.pybuiltins.list) ? [pyconvert(AbstractProgramResult, rt) for rt in x_jaqcd.results] : AbstractProgramResult[]
    bris         = pyisinstance(x_jaqcd.basis_rotation_instructions, PythonCall.pybuiltins.list) ? [pyconvert(Instruction, ix) for ix in x_jaqcd.basis_rotation_instructions] : Instruction[]
    prog         = Braket.Program(Braket.header_dict[Braket.Program], instructions, results, bris)
    PythonCall.pyconvert_return(prog)
end
