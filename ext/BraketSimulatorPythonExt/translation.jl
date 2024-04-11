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
for (rt, brt, fn) in ((:(Braket.IR.Sample), :(Braket.Sample), :jl_convert_sample),
                      (:(Braket.IR.Expectation), :(Braket.Expectation), :jl_convert_expectation),
                      (:(Braket.IR.Variance), :(Braket.Variance), :jl_convert_variance),
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
        function $fn(::Type{Braket.Result}, x::Py)
            jl_obs = Braket.StructTypes.constructfrom(Braket.Observables.Observable, convert_ir_obs(x.observable))
            jl_ts  = convert_ir_target(x.targets)
            jl_rt  = $brt(jl_obs, jl_ts)
            PythonCall.pyconvert_return(jl_rt)
        end
        function $fn(::Type{$brt}, x::Py)
            jl_obs = Braket.StructTypes.constructfrom(Braket.Observables.Observable, convert_ir_obs(x.observable))
            jl_ts  = convert_ir_target(x.targets)
            jl_rt  = $brt(jl_obs, jl_ts)
            PythonCall.pyconvert_return(jl_rt)
        end
    end
end

for (rt, brt, fn) in ((:(Braket.IR.Probability), :(Braket.Probability), :jl_convert_probability),
                      (:(Braket.IR.DensityMatrix), :(Braket.DensityMatrix), :jl_convert_densitymatrix),
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
        function $fn(::Type{Braket.Result}, x::Py)
            jl_ts  = convert_ir_target(x.targets)
            PythonCall.pyconvert_return($brt(jl_ts))
        end
        function $fn(::Type{$brt}, x::Py)
            jl_ts  = convert_ir_target(x.targets)
            PythonCall.pyconvert_return($brt(jl_ts))
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

function jl_convert_amplitude(::Type{Braket.Amplitude}, x::Py)
    jl_rt = Braket.Amplitude([pyconvert(String, s) for s in x.states])
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_amplitude(::Type{Braket.Result}, x::Py)
    jl_rt = Braket.Amplitude([pyconvert(String, s) for s in x.states])
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
function jl_convert_statevector(::Type{Braket.StateVector}, x::Py)
    jl_rt = Braket.StateVector()
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_statevector(::Type{Braket.Result}, x::Py)
    jl_rt = Braket.StateVector()
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

for (fn, jl_typ) in ((:jl_convert_sim_identity, :(Braket.I)),
                     (:jl_convert_sim_hadamard, :H),
                     (:jl_convert_sim_paulix, :X),
                     (:jl_convert_sim_pauliy, :Y),
                     (:jl_convert_sim_pauliz, :Z),
                     (:jl_convert_sim_cv, :CV),
                     (:jl_convert_sim_cx, :CNot),
                     (:jl_convert_sim_cy, :CY),
                     (:jl_convert_sim_cz, :CZ),
                     (:jl_convert_sim_v, :V),
                     (:jl_convert_sim_vi, :Vi),
                     (:jl_convert_sim_t, :T),
                     (:jl_convert_sim_ti, :Ti),
                     (:jl_convert_sim_s, :S),
                     (:jl_convert_sim_si, :Si),
                     (:jl_convert_sim_ecr, :ECR),
                     (:jl_convert_sim_swap, :Swap),
                     (:jl_convert_sim_iswap, :ISwap),
                     (:jl_convert_sim_cswap, :CSwap),
                     (:jl_convert_sim_ccnot, :CCNot),
                    )
    @eval begin
        function $fn(::Type{Braket.Instruction}, x::Py)
            jl_x       = $jl_typ()
            jl_targets = Int[pyconvert(Int, t) for t in x.targets]
            exponent   = pyconvert(Any, x._power)
            jl_x       = !isone(exponent) ? jl_x ^ exponent : jl_x
            jl_op      = Control(jl_x, tuple((pyconvert(Int, c) for c in x._ctrl_modifiers)...))
            PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
        end
    end
end
for (fn, jl_typ) in ((:jl_convert_sim_identity, :(Braket.Observables.I)),
                     (:jl_convert_sim_hadamard, :(Braket.Observables.H)),
                     (:jl_convert_sim_paulix, :(Braket.Observables.X)),
                     (:jl_convert_sim_pauliy, :(Braket.Observables.Y)),
                     (:jl_convert_sim_pauliz, :(Braket.Observables.Z)),
                    )
    @eval begin
        $fn(::Type{$jl_typ}, x::Py) = PythonCall.pyconvert_return($jl_typ())
        $fn(::Type{Braket.Observables.Observable}, x::Py) = PythonCall.pyconvert_return($jl_typ())
    end
end
jl_convert_sim_tensorproduct(::Type{Braket.Observables.TensorProduct}, x::Py) = PythonCall.pyconvert_return(Braket.Observables.TensorProduct([pyconvert(Braket.Observables.Observable, f) for f in x._factors]))
jl_convert_sim_tensorproduct(::Type{Braket.Observables.Observable}, x::Py)    = PythonCall.pyconvert_return(Braket.Observables.TensorProduct([pyconvert(Braket.Observables.Observable, f) for f in x._factors]))

jl_convert_sim_hermitian(::Type{Braket.Observables.HermitianObservable}, x::Py) = PythonCall.pyconvert_return(Braket.Observables.HermitianObservable(pyconvert(Matrix{ComplexF64}, x._matrix)))
jl_convert_sim_hermitian(::Type{Braket.Observables.Observable}, x::Py)    = PythonCall.pyconvert_return(Braket.Observables.HermitianObservable(pyconvert(Matrix{ComplexF64}, x._matrix)))

for (fn, jl_typ) in ((:jl_convert_sim_gpi, :GPi),
                     (:jl_convert_sim_gpi2, :GPi2),
                     (:jl_convert_sim_phaseshift, :PhaseShift),
                     (:jl_convert_sim_cphaseshift, :CPhaseShift),
                     (:jl_convert_sim_cphaseshift00, :CPhaseShift00),
                     (:jl_convert_sim_cphaseshift01, :CPhaseShift01),
                     (:jl_convert_sim_cphaseshift10, :CPhaseShift10),
                     (:jl_convert_sim_pswap, :PSwap),
                     (:jl_convert_sim_rx, :Rx),
                     (:jl_convert_sim_ry, :Ry),
                     (:jl_convert_sim_rz, :Rz),
                     (:jl_convert_sim_xx, :XX),
                     (:jl_convert_sim_xy, :XY),
                     (:jl_convert_sim_yy, :YY),
                     (:jl_convert_sim_zz, :ZZ),
                    )
    @eval begin
        function $fn(::Type{Braket.Instruction}, x::Py)
            jl_x       = $jl_typ(pyconvert(Union{Float64, FreeParameter}, x._angle))
            jl_targets = Int[pyconvert(Int, t) for t in x.targets]
            exponent   = pyconvert(Any, x._power)
            jl_x       = !isone(exponent) ? jl_x ^ exponent : jl_x
            jl_op      = Control(jl_x, tuple((pyconvert(Int, c) for c in x._ctrl_modifiers)...))
            PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
        end
    end
end
for (fn, jl_typ) in ((:jl_convert_sim_bitflip, :BitFlip),
                     (:jl_convert_sim_phaseflip, :PhaseFlip),
                     (:jl_convert_sim_depolarizing, :Depolarizing),
                     (:jl_convert_sim_twoqubitdepolarizing, :TwoQubitDepolarizing),
                     (:jl_convert_sim_twoqubitdephasing, :TwoQubitDephasing),
                    )
    @eval begin
        function $fn(::Type{Braket.Instruction}, x::Py)
            jl_op      = $jl_typ(pyconvert(Union{Float64, FreeParameter}, x._probability))
            jl_targets = Int[pyconvert(Int, t) for t in x.targets]
            PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
        end
    end
end

for (fn, jl_typ) in ((:jl_convert_sim_amplitudedamping, :AmplitudeDamping),
                     (:jl_convert_sim_phasedamping, :PhaseDamping),
                    )
    @eval begin
        function $fn(::Type{Braket.Instruction}, x::Py)
            jl_op      = $jl_typ(pyconvert(Union{Float64, FreeParameter}, x._gamma))
            jl_targets = Int[pyconvert(Int, t) for t in x.targets]
            PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
        end
    end
end
function jl_convert_sim_paulichannel(::Type{Braket.Instruction}, x::Py)
    jl_op      = PauliChannel(pyconvert(Union{Float64, FreeParameter}, x._probX), pyconvert(Union{Float64, FreeParameter}, x._probY), pyconvert(Union{Float64, FreeParameter}, x._probZ))
    jl_targets = Int[pyconvert(Int, t) for t in x.targets]
    PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
end
function jl_convert_sim_generalizedamplitudedamping(::Type{Braket.Instruction}, x::Py)
    jl_op      = GeneralizedAmplitudeDamping(pyconvert(Union{Float64, FreeParameter}, x._probability), pyconvert(Union{Float64, FreeParameter}, x._gamma))
    jl_targets = Int[pyconvert(Int, t) for t in x.targets]
    PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
end
for (fn, jl_typ, a1, a2, a3) in ((:jl_convert_sim_ms, :MS, "_angle_1", "_angle_2", "_angle_3"),
                                 (:jl_convert_sim_u, :U, "_theta", "_phi", "_lambda"),
                                )
    @eval begin
        function $fn(::Type{Braket.Instruction}, x::Py)
            angle      = (pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a1))),
                          pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a2))),
                          pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a3))))
            jl_x       = $jl_typ(angle)
            jl_targets = Int[pyconvert(Int, t) for t in x.targets]
            exponent   = pyconvert(Any, x._power)
            jl_x       = !isone(exponent) ? jl_x ^ exponent : jl_x
            jl_op      = Control(jl_x, tuple((pyconvert(Int, c) for c in x._ctrl_modifiers)...))
            PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
        end
    end
end

jl_convert_sympy_Pi(::Type{Float64}, x::Py) = PythonCall.pyconvert_return(convert(Float64, π))
jl_convert_sympy_E(::Type{Float64}, x::Py) = PythonCall.pyconvert_return(convert(Float64, ℯ))

function jl_convert_sympy_Mul(::Type{Float64}, x::Py)
    jl_args = [pyconvert(Float64, a) for a in x.args]
    val     = prod(jl_args, init=1.0)
    PythonCall.pyconvert_return(val)
end

function jl_convert_sim_gphase(::Type{Braket.Instruction}, x::Py)
    angle      = pyconvert(Union{FreeParameter, Float64}, x._angle)
    jl_targets = Int[pyconvert(Int, t) for t in x.targets]
    jl_op       = MultiQubitPhaseShift{length(jl_targets)}(angle)
    PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
end
function jl_convert_sim_unitary(::Type{Braket.Instruction}, x::Py)
    jl_mat     = pyconvert(Matrix{ComplexF64}, x._matrix)
    jl_x       = Unitary(jl_mat)
    jl_targets = Int[pyconvert(Int, t) for t in x.targets]
    exponent   = pyconvert(Any, x._power)
    jl_x       = !isone(exponent) ? jl_x ^ exponent : jl_x
    jl_op      = Control(jl_x, tuple((pyconvert(Int, c) for c in x._ctrl_modifiers)...))
    PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
end
function jl_convert_sim_kraus(::Type{Braket.Instruction}, x::Py)
    jl_mats    = [pyconvert(Matrix{ComplexF64}, m) for m in x._matrices]
    jl_op      = Kraus(jl_mats)
    jl_targets = Int[pyconvert(Int, t) for t in x.targets]
    PythonCall.pyconvert_return(Instruction(jl_op, jl_targets))
end

for (rt, fn) in ((:(Braket.Sample), :jl_convert_sim_sample),
                 (:(Braket.Expectation), :jl_convert_sim_expectation),
                 (:(Braket.Variance), :jl_convert_sim_variance),
                )
    @eval begin
        function $fn(::Type{$rt}, x::Py)
            jl_obs = pyconvert(Braket.Observables.Observable, x._observable)
            jl_ts  = convert_ir_target(x._targets)
            jl_rt  = $rt(jl_obs, jl_ts)
            PythonCall.pyconvert_return(jl_rt)
        end
        function $fn(::Type{Braket.Result}, x::Py)
            jl_obs = pyconvert(Braket.Observables.Observable, x._observable)
            jl_ts  = convert_ir_target(x._targets)
            jl_rt  = $rt(jl_obs, jl_ts)
            PythonCall.pyconvert_return(jl_rt)
        end
    end
end

for (rt, fn) in ((:(Braket.Probability), :jl_convert_sim_probability),
                 (:(Braket.DensityMatrix), :jl_convert_sim_densitymatrix),
                )
    @eval begin
        $fn(::Type{$rt}, x::Py) = PythonCall.pyconvert_return($rt(convert_ir_target(x.targets)))
        $fn(::Type{Braket.Result}, x::Py) = PythonCall.pyconvert_return($rt(convert_ir_target(x.targets)))
    end
end

function jl_convert_sim_amplitude(::Type{Braket.Amplitude}, x::Py)
    jl_rt = Braket.Amplitude([pyconvert(String, s) for s in x._states])
    PythonCall.pyconvert_return(jl_rt)
end
function jl_convert_sim_amplitude(::Type{Braket.Result}, x::Py)
    jl_rt = Braket.Amplitude([pyconvert(String, s) for s in x._states])
    PythonCall.pyconvert_return(jl_rt)
end

function jl_convert_sim_adjointgradient(::Type{Braket.AdjointGradient}, x::Py)
    error("not implemented yet!")
end
function jl_convert_sim_adjointgradient(::Type{Braket.Result}, x::Py)
    error("not implemented yet!")
end

jl_convert_sim_statevector(::Type{Braket.StateVector}, x::Py) = PythonCall.pyconvert_return(Braket.StateVector())
jl_convert_sim_statevector(::Type{Braket.Result}, x::Py)      = PythonCall.pyconvert_return(Braket.StateVector())

function jl_convert_sim_circuit(::Type{Braket.Circuit}, x)
    instructions = [pyconvert(Instruction, ix) for ix in x.instructions]
    results      = [pyconvert(Result, rt) for rt in x.results]
    prog         = Braket.Circuit()
    for ix in instructions
        Braket.add_instruction!(prog, ix)
    end
    for rt in results 
        push!(prog.result_types, rt)
    end
    PythonCall.pyconvert_return(prog)
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
