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

convert_targets(ts::Py)          = Int[pyconvert(Int, t) for t in ts]
convert_targets(::Nothing)       = PythonCall.pybuiltins.None
convert_targets(ts::Vector{Int}) = pylist(ts)
function convert_ir_target(t)
    pyis(t, pybuiltins.None) && return nothing
    return convert_targets(t)
end

for (rt, brt, fn) in ((:(Braket.IR.Sample), :(Braket.Sample), :jl_convert_sample),
                      (:(Braket.IR.Expectation), :(Braket.Expectation), :jl_convert_expectation),
                      (:(Braket.IR.Variance), :(Braket.Variance), :jl_convert_variance),
                     )
    @eval begin
        function $fn(t::Type{$rt}, x::Py)
            jl_obs = convert_ir_obs(x.observable)
            jl_ts  = convert_ir_target(x.targets)
            jl_rt  = t(jl_obs, jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
        $fn(::Type{AbstractProgramResult}, x::Py) = $fn($rt, x)
        function $fn(t::Type{$brt}, x::Py)
            jl_obs = Braket.StructTypes.constructfrom(Braket.Observables.Observable, convert_ir_obs(x.observable))
            PythonCall.pyconvert_return(t(jl_obs, convert_ir_target(x.targets)))
        end
        $fn(::Type{Braket.Result}, x::Py) = $fn($brt, x::Py)
    end
end

for (rt, brt, fn) in ((:(Braket.IR.Probability), :(Braket.Probability), :jl_convert_probability),
                      (:(Braket.IR.DensityMatrix), :(Braket.DensityMatrix), :jl_convert_densitymatrix),
                     )
    @eval begin
        function $fn(t::Type{$rt}, x::Py)
            jl_ts = convert_ir_target(x.targets)
            jl_rt = t(jl_ts, pyconvert(String, x.type))
            PythonCall.pyconvert_return(jl_rt)
        end
        $fn(::Type{AbstractProgramResult}, x::Py) = $fn($rt, x)
        $fn(::Type{$brt}, x::Py) = PythonCall.pyconvert_return($brt(convert_ir_target(x.targets)))
        $fn(::Type{Braket.Result}, x::Py) = $fn($brt, x::Py)
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
function jl_convert_oqprogram(t::Type{OpenQasmProgram}, x::Py)
    bsh    = pyconvert(Braket.braketSchemaHeader, x.braketSchemaHeader)
    source = pyconvert(String, x.source)
    inputs = inputs_to_jl(x.inputs) 
    PythonCall.pyconvert_return(t(bsh, source, inputs))
end

function jl_convert_amplitude(t::Type{Braket.IR.Amplitude}, x::Py)
    states = map(state->pyconvert(String, state), x.states)
    jl_rt  = t(states, pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
jl_convert_amplitude(::Type{AbstractProgramResult}, x::Py) = jl_convert_amplitude(Braket.IR.Amplitude, x::Py)

function jl_convert_amplitude(t::Type{Braket.Amplitude}, x::Py)
    states = map(state->pyconvert(String, state), x.states)
    PythonCall.pyconvert_return(t(states))
end
jl_convert_amplitude(::Type{Braket.Result}, x::Py) = jl_convert_amplitude(Braket.Amplitude, x::Py)

function jl_convert_statevector(t::Type{Braket.IR.StateVector}, x::Py)
    jl_rt = t(pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
jl_convert_statevector(::Type{AbstractProgramResult}, x::Py) = jl_convert_statevector(Braket.IR.StateVector, x::Py)
jl_convert_statevector(::Type{Braket.StateVector}, x::Py) = PythonCall.pyconvert_return(Braket.StateVector())
jl_convert_statevector(::Type{Braket.Result}, x::Py) = jl_convert_statevector(Braket.StateVector, x::Py)

# exclude adjoint gradient translation from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
function jl_convert_adjointgradient(t::Type{Braket.IR.AdjointGradient}, x::Py)
    jl_targets = pyis(x.targets, pybuiltins.None) ? nothing : [convert_targets(t_) for t_ in t]
    jl_params = pyis(x.targets, pybuiltins.None) ? nothing : [pyconvert(String, p) for p in x.parameters]
    jl_obs = pyisinstance(x.observable, pybuiltins.str) ? pyconvert(String, x.observable) : [convert_ir_obs(o) for o in x.observable] 
    jl_rt = t(jl_params, jl_obs, jl_targets, pyconvert(String, x.type))
    PythonCall.pyconvert_return(jl_rt)
end
jl_convert_adjointgradient(::Type{AbstractProgramResult}, x::Py) = jl_convert_adjointgradient(Braket.IR.AdjointGradient, x::Py)
# COV_EXCL_STOP

jl_convert_kraus(::Type{Kraus}, x::Py)::Kraus = PythonCall.pyconvert_return(Kraus([convert_ir_matrix(m) for m in x.matrices]))
jl_convert_kraus(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Kraus([convert_ir_matrix(m) for m in x.matrices]), convert_targets(x.targets)))

jl_convert_unitary(::Type{Unitary}, x::Py)::Unitary = PythonCall.pyconvert_return(Unitary(convert_ir_matrix(x.matrix)))
jl_convert_unitary(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Unitary(convert_ir_matrix(x.matrix)), convert_targets(x.targets)))

jl_convert_cswap(::Type{CSwap}, x::Py)::CSwap = PythonCall.pyconvert_return(CSwap())
jl_convert_cswap(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(CSwap(), vcat(pyconvert(Int, x.control), convert_targets(x.targets))))
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
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{Float64, FreeParameter}, x.probability)), convert_targets(x.targets)))
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
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{Float64, FreeParameter}, x.angle)))
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
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Union{FreeParameter, Float64}, x.angle)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Union{FreeParameter, Float64}, x.angle)), convert_targets(x.targets)))
    end
end

for (g, fn) in ((:StartVerbatimBox, :jl_convert_startverbatim),
                (:EndVerbatimBox, :jl_convert_endverbatim),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g())
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g()))
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

function jl_convert_ix(t::Type{Instruction}, x::Py)::Instruction
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
    sts     = merge(Braket.StructTypes.subtypes(Gate), Braket.StructTypes.subtypes(Noise), Braket.StructTypes.subtypes(CompilerDirective))
    py_name = pyconvert(String, pygetattr(x, "type"))
    jl_type = sts[Symbol(py_name)]
    # bad, can we come up with something more generic
    ix = t(PythonCall.pyconvert(jl_type, x), full_targets)
    PythonCall.pyconvert_return(ix)
end

function jl_convert_bsh(t::Type{Braket.braketSchemaHeader}, x::Py)
    name    = pyconvert(String, x.name)
    version = pyconvert(String, x.version)
    PythonCall.pyconvert_return(t(name, version))
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

function _handle_modifiers(raw_operation, x::Py)
    exponent        = pyconvert(Any, x._power)
    power_operation = !isone(exponent) ? raw_operation ^ exponent : raw_operation
    return Control(power_operation, tuple((pyconvert(Int, c) for c in x._ctrl_modifiers)...))
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
        function $fn(t::Type{Braket.Instruction}, x::Py)
            jl_operation = _handle_modifiers($jl_typ(), x)
            PythonCall.pyconvert_return(t(jl_operation, convert_targets(x.targets)))
        end
    end
end
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
        function $fn(t::Type{Braket.Instruction}, x::Py)
            jl_operation = _handle_modifiers($jl_typ(pyconvert(Union{Float64, FreeParameter}, x._angle)), x)
            PythonCall.pyconvert_return(t(jl_operation, convert_targets(x.targets)))
        end
    end
end
function jl_convert_sim_prx(t::Type{Braket.Instruction}, x::Py)
    angle = (pyconvert(Union{Float64, FreeParameter}, getproperty(x, :_angle_1)),
             pyconvert(Union{Float64, FreeParameter}, getproperty(x, :_angle_2)),
            )
    jl_operation = _handle_modifiers(PRx(angle), x)
    PythonCall.pyconvert_return(t(jl_operation, convert_targets(x.targets)))
end
for (fn, jl_typ, a1, a2, a3) in ((:jl_convert_sim_ms, :MS, "_angle_1", "_angle_2", "_angle_3"),
                                 (:jl_convert_sim_u, :U, "_theta", "_phi", "_lambda"),
                                )
    @eval begin
        function $fn(t::Type{Braket.Instruction}, x::Py)
            angle      = (pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a1))),
                          pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a2))),
                          pyconvert(Union{Float64, FreeParameter}, getproperty(x, Symbol($a3))))
            jl_operation = _handle_modifiers($jl_typ(angle), x)
            PythonCall.pyconvert_return(t(jl_operation, convert_targets(x.targets)))
        end
    end
end
function jl_convert_sim_unitary(t::Type{Braket.Instruction}, x::Py)
    jl_matrix    = pyconvert(Matrix{ComplexF64}, x._matrix)
    jl_operation = _handle_modifiers(Unitary(jl_matrix), x)
    PythonCall.pyconvert_return(t(jl_operation, convert_targets(x.targets)))
end
for (fn, jl_typ) in ((:jl_convert_sim_bitflip, :BitFlip),
                     (:jl_convert_sim_phaseflip, :PhaseFlip),
                     (:jl_convert_sim_depolarizing, :Depolarizing),
                     (:jl_convert_sim_twoqubitdepolarizing, :TwoQubitDepolarizing),
                     (:jl_convert_sim_twoqubitdephasing, :TwoQubitDephasing),
                    )
    @eval begin
        function $fn(t::Type{Braket.Instruction}, x::Py)
            jl_op = $jl_typ(pyconvert(Union{Float64, FreeParameter}, x._probability))
            PythonCall.pyconvert_return(t(jl_op, convert_targets(x.targets)))
        end
    end
end

for (fn, jl_typ) in ((:jl_convert_sim_amplitudedamping, :AmplitudeDamping),
                     (:jl_convert_sim_phasedamping, :PhaseDamping),
                    )
    @eval begin
        function $fn(t::Type{Braket.Instruction}, x::Py)
            jl_op = $jl_typ(pyconvert(Union{Float64, FreeParameter}, x._gamma))
            PythonCall.pyconvert_return(t(jl_op, convert_targets(x.targets)))
        end
    end
end
function jl_convert_sim_paulichannel(t::Type{Braket.Instruction}, x::Py)
    jl_op = PauliChannel(pyconvert(Union{Float64, FreeParameter}, x._probX),
                         pyconvert(Union{Float64, FreeParameter}, x._probY),
                         pyconvert(Union{Float64, FreeParameter}, x._probZ)
                        )
    PythonCall.pyconvert_return(t(jl_op, convert_targets(x.targets)))
end
function jl_convert_sim_generalizedamplitudedamping(t::Type{Braket.Instruction}, x::Py)
    jl_op = GeneralizedAmplitudeDamping(pyconvert(Union{Float64, FreeParameter}, x._probability),
                                        pyconvert(Union{Float64, FreeParameter}, x._gamma),
                                       )
    PythonCall.pyconvert_return(t(jl_op, convert_targets(x.targets)))
end

jl_convert_sympy_Pi(::Type{Float64}, x::Py) = PythonCall.pyconvert_return(convert(Float64, π))
jl_convert_sympy_E(::Type{Float64}, x::Py) = PythonCall.pyconvert_return(convert(Float64, ℯ))

jl_convert_sympy_Mul(::Type{Float64}, x::Py) = PythonCall.pyconvert_return(mapreduce(arg->pyconvert(Float64, arg), *, x.args, init=1.0))

function jl_convert_sim_gphase(t::Type{Braket.Instruction}, x::Py)
    angle      = pyconvert(Union{FreeParameter, Float64}, x._angle)
    jl_targets = convert_targets(x.targets)
    jl_op      = MultiQubitPhaseShift{length(jl_targets)}(angle)
    PythonCall.pyconvert_return(t(jl_op, jl_targets))
end

function jl_convert_sim_kraus(t::Type{Braket.Instruction}, x::Py)
    jl_op      = Kraus([pyconvert(Matrix{ComplexF64}, m) for m in x._matrices])
    PythonCall.pyconvert_return(t(jl_op, convert_targets(x.targets)))
end

for (rt, fn) in ((:(Braket.Sample), :jl_convert_sim_sample),
                 (:(Braket.Expectation), :jl_convert_sim_expectation),
                 (:(Braket.Variance), :jl_convert_sim_variance),
                )
    @eval begin
        function $fn(t::Type{$rt}, x::Py)
            jl_obs = pyconvert(Braket.Observables.Observable, x._observable)
            PythonCall.pyconvert_return(t(jl_obs, convert_ir_target(x._targets)))
        end
    end
end

for (rt, fn) in ((:(Braket.Probability), :jl_convert_sim_probability),
                 (:(Braket.DensityMatrix), :jl_convert_sim_densitymatrix),
                )
    @eval begin
        $fn(::Type{$rt}, x::Py) = PythonCall.pyconvert_return($rt(convert_ir_target(x.targets)))
    end
end

function jl_convert_sim_amplitude(t::Type{Braket.Amplitude}, x::Py)
    states = map(state->pyconvert(String, state), x._states)
    PythonCall.pyconvert_return(t(states))
end

# exclude adjoint gradient translation from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
jl_convert_sim_adjointgradient(::Type{Braket.AdjointGradient}, x::Py) = error("not implemented yet!")
jl_convert_sim_adjointgradient(::Type{Braket.Result}, x::Py) = error("not implemented yet!")
# COV_EXCL_STOP
jl_convert_sim_statevector(::Type{Braket.StateVector}, x::Py) = PythonCall.pyconvert_return(Braket.StateVector())

function jl_convert_sim_circuit(t::Type{Braket.Circuit}, x)
    instructions = map(ix->pyconvert(Instruction, ix), x.instructions)
    results      = map(rt->pyconvert(Result, rt), x.results)
    prog         = t()
    foreach(ix->Braket.add_instruction!(prog, ix), instructions)
    foreach(rt->push!(prog.result_types, rt), results)
    PythonCall.pyconvert_return(prog)
end

function jl_convert_program(t::Type{Braket.IR.Program}, x_jaqcd)
    instructions = map(ix->pyconvert(Instruction, ix), x_jaqcd.instructions)
    results      = pyisinstance(x_jaqcd.results, PythonCall.pybuiltins.list) ? [pyconvert(AbstractProgramResult, rt) for rt in x_jaqcd.results] : AbstractProgramResult[]
    bris         = pyisinstance(x_jaqcd.basis_rotation_instructions, PythonCall.pybuiltins.list) ? [pyconvert(Instruction, ix) for ix in x_jaqcd.basis_rotation_instructions] : Instruction[]
    prog         = t(Braket.header_dict[Braket.Program], instructions, results, bris)
    PythonCall.pyconvert_return(prog)
end
jl_convert_circuit(::Type{Braket.IR.Program}, x) = jl_convert_program(Braket.IR.Program, x._to_jaqcd())

py_obs(o::String) = pylist([pystr(o)])
function py_obs(obs::Vector)
    raw_obs = map(obs) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    return pylist(raw_obs)
end
Py(r::Braket.IR.Sample)        = braket[].ir.jaqcd.results.Sample(targets=convert_targets(r.targets), observable=py_obs(r.observable), type=pystr("sample"))
Py(r::Braket.IR.Expectation)   = braket[].ir.jaqcd.results.Expectation(targets=convert_targets(r.targets), observable=py_obs(r.observable), type=pystr("expectation"))
Py(r::Braket.IR.Variance)      = braket[].ir.jaqcd.results.Variance(observable=py_obs(r.observable), targets=convert_targets(r.targets), type=pystr("variance"))
Py(r::Braket.IR.Amplitude)     = braket[].ir.jaqcd.results.Amplitude(states=pylist(pystr(s) for s in r.states), type=pystr("amplitude"))
Py(r::Braket.IR.StateVector)   = braket[].ir.jaqcd.results.StateVector(type=pystr("statevector"))
Py(r::Braket.IR.DensityMatrix) = braket[].ir.jaqcd.results.DensityMatrix(targets=convert_targets(r.targets), type=pystr("densitymatrix"))
Py(r::Braket.IR.Probability)   = braket[].ir.jaqcd.results.Probability(targets=convert_targets(r.targets), type=pystr("probability"))
# exclude adjoint gradient translation from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
Py(r::Braket.IR.AdjointGradient) =  braket[].ir.jaqcd.results.AdjointGradient(targets=convert_targets(r.targets), observable=pylist(r.observable), parameters=pylist(pystr(p) for p in r.parameters), type=pystr("adjoint_gradient"))
# COV_EXCL_STOP

function Py(rt::Braket.ResultTypeValue)
    py_typ = Py(rt.type)
    py_val = if rt.value isa Dict
        pydict(rt.value)
    elseif rt.value isa Float64
        rt.value
    elseif rt.value isa Vector{Vector{Vector{Float64}}}
        pylist(pylist(pycomplex(v_...) for v_ in v) for v in rt.value)
    else
        pylist(rt.value)
    end
    return braket[].task_result.gate_model_task_result_v1.ResultTypeValue(type=py_typ, value=py_val)
end

Py(op::Braket.H, targets) = braket[].ir.jaqcd.instructions.H(target=targets[1], type=pystr("h"))
Py(op::Braket.X, targets) = braket[].ir.jaqcd.instructions.X(target=targets[1], type=pystr("x"))
Py(op::Braket.Y, targets) = braket[].ir.jaqcd.instructions.Y(target=targets[1], type=pystr("y"))
Py(op::Braket.Z, targets) = braket[].ir.jaqcd.instructions.Z(target=targets[1], type=pystr("z"))
Py(op::Braket.I, targets) = braket[].ir.jaqcd.instructions.I(target=targets[1], type=pystr("i"))
Py(op::Braket.T, targets) = braket[].ir.jaqcd.instructions.T(target=targets[1], type=pystr("t"))
Py(op::Braket.V, targets) = braket[].ir.jaqcd.instructions.V(target=targets[1], type=pystr("v"))
Py(op::Braket.S, targets) = braket[].ir.jaqcd.instructions.S(target=targets[1], type=pystr("s"))
Py(op::Braket.Ti, targets) = braket[].ir.jaqcd.instructions.Ti(target=targets[1], type=pystr("ti"))
Py(op::Braket.Si, targets) = braket[].ir.jaqcd.instructions.Si(target=targets[1], type=pystr("si"))
Py(op::Braket.Vi, targets) = braket[].ir.jaqcd.instructions.Vi(target=targets[1], type=pystr("vi"))
Py(op::Braket.Rx, targets) = braket[].ir.jaqcd.instructions.Rx(target=targets[1], angle=op.angle[1], type=pystr("rx"))
Py(op::Braket.Ry, targets) = braket[].ir.jaqcd.instructions.Ry(target=targets[1], angle=op.angle[1], type=pystr("ry"))
Py(op::Braket.Rz, targets) = braket[].ir.jaqcd.instructions.Rz(target=targets[1], angle=op.angle[1], type=pystr("rz"))
Py(op::Braket.XX, targets) = braket[].ir.jaqcd.instructions.XX(targets=pylist(targets), angle=op.angle[1], type=pystr("xx"))
Py(op::Braket.YY, targets) = braket[].ir.jaqcd.instructions.YY(targets=pylist(targets), angle=op.angle[1], type=pystr("yy"))
Py(op::Braket.ZZ, targets) = braket[].ir.jaqcd.instructions.ZZ(targets=pylist(targets), angle=op.angle[1], type=pystr("zz"))
Py(op::Braket.XY, targets) = braket[].ir.jaqcd.instructions.XY(targets=pylist(targets), angle=op.angle[1], type=pystr("xy"))
Py(op::Braket.Swap, targets) = braket[].ir.jaqcd.instructions.Swap(targets=pylist(targets), type=pystr("swap"))
Py(op::Braket.ISwap, targets) = braket[].ir.jaqcd.instructions.ISwap(targets=pylist(targets), type=pystr("iswap"))
Py(op::Braket.PSwap, targets) = braket[].ir.jaqcd.instructions.PSwap(targets=pylist(targets), angle=op.angle[1], type=pystr("pswap"))
Py(op::Braket.CSwap, targets) = braket[].ir.jaqcd.instructions.CSwap(control=targets[1], targets=pylist(targets[2:end]), type=pystr("cswap"))
Py(op::Braket.CZ, targets) = braket[].ir.jaqcd.instructions.CZ(control=targets[1], target=targets[2], type=pystr("cz"))
Py(op::Braket.CY, targets) = braket[].ir.jaqcd.instructions.CY(control=targets[1], target=targets[2], type=pystr("cy"))
Py(op::Braket.CV, targets) = braket[].ir.jaqcd.instructions.CV(control=targets[1], target=targets[2], type=pystr("cv"))
Py(op::Braket.CNot, targets) = braket[].ir.jaqcd.instructions.CNot(control=targets[1], target=targets[2], type=pystr("cnot"))
Py(op::Braket.CCNot, targets) = braket[].ir.jaqcd.instructions.CCNot(controls=pylist(targets[1:2]), target=targets[3], type=pystr("ccnot"))
Py(op::Braket.PhaseShift, targets) = braket[].ir.jaqcd.instructions.PhaseShift(target=targets[1], angle=op.angle[1], type=pystr("phaseshift"))
Py(op::Braket.CPhaseShift, targets) = braket[].ir.jaqcd.instructions.CPhaseShift(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift"))
Py(op::Braket.CPhaseShift00, targets) = braket[].ir.jaqcd.instructions.CPhaseShift00(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift00"))
Py(op::Braket.CPhaseShift01, targets) = braket[].ir.jaqcd.instructions.CPhaseShift01(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift01"))
Py(op::Braket.CPhaseShift10, targets) = braket[].ir.jaqcd.instructions.CPhaseShift10(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift10"))
Py(op::Braket.ECR, targets) = braket[].ir.jaqcd.instructions.ECR(targets=pylist(targets), type=pystr("ecr"))
Py(op::Braket.BitFlip, targets) = braket[].ir.jaqcd.instructions.BitFlip(target=targets[1], probability=op.probability, type=pystr("bit_flip"))
Py(op::Braket.PhaseFlip, targets) = braket[].ir.jaqcd.instructions.PhaseFlip(target=targets[1], probability=op.probability, type=pystr("phase_flip"))
Py(op::Braket.PauliChannel, targets) = braket[].ir.jaqcd.instructions.PauliChannel(target=targets[1], probX=op.probX, probY=op.probY, probZ=op.probZ, type=pystr("pauli_channel"))
Py(op::Braket.MultiQubitPauliChannel, targets) = braket[].ir.jaqcd.instructions.MultiQubitPauliChannel(target=pylist(targets), probabilities=pydict(op.probabilities), type=pystr("multi_qubit_pauli_channel"))
Py(op::Braket.Depolarizing, targets) = braket[].ir.jaqcd.instructions.Depolarizing(target=targets[1], probability=op.probability, type=pystr("depolarizing"))
Py(op::Braket.TwoQubitDepolarizing, targets) = braket[].ir.jaqcd.instructions.TwoQubitDepolarizing(targets=pylist(targets), probability=op.probability, type=pystr("two_qubit_depolarizing"))
Py(op::Braket.TwoQubitDephasing, targets) = braket[].ir.jaqcd.instructions.TwoQubitDephasing(targets=pylist(targets), probability=op.probability, type=pystr("two_qubit_dephasing"))
Py(op::Braket.AmplitudeDamping, targets) = braket[].ir.jaqcd.instructions.AmplitudeDamping(target=targets[1], gamma=op.gamma, type=pystr("amplitude_damping"))
Py(op::Braket.PhaseDamping, targets) = braket[].ir.jaqcd.instructions.PhaseDamping(target=targets[1], gamma=op.gamma, type=pystr("phase_damping"))
Py(op::Braket.GeneralizedAmplitudeDamping, targets) = braket[].ir.jaqcd.instructions.GeneralizedAmplitudeDamping(target=targets[1], gamma=op.gamma, probability=op.probability, type=pystr("generalized_amplitude_damping"))
Py(op::Braket.StartVerbatimBox) = braket[].ir.jaqcd.instructions.StartVerbatimBox(type=pystr("start_verbatim_box"))
Py(op::Braket.EndVerbatimBox) = braket[].ir.jaqcd.instructions.EndVerbatimBox(type=pystr("end_verbatim_box"))
Py(op::Braket.StartVerbatimBox, targets) = braket[].ir.jaqcd.instructions.StartVerbatimBox(type=pystr("start_verbatim_box"))
Py(op::Braket.EndVerbatimBox, targets) = braket[].ir.jaqcd.instructions.EndVerbatimBox(type=pystr("end_verbatim_box"))
function Py(op::Braket.Unitary, targets)
    raw_mat = Braket.complex_matrix_to_ir(op.matrix)
    py_mat = pylist(pylist(pylist(v__) for v__ in v_) for v_ in raw_mat)
    return braket[].ir.jaqcd.instructions.Unitary(targets=pylist(targets), matrix=py_mat, type=pystr("unitary"))
end
function Py(op::Braket.Kraus, targets)
    raw_mats = map(Braket.complex_matrix_to_ir, op.matrices)
    py_mat = pylist(pylist(pylist(pylist(v__) for v__ in v_) for v_ in raw_mat) for raw_mat in raw_mats)
    return braket[].ir.jaqcd.instructions.Kraus(targets=pylist(targets), matrices=py_mat, type=pystr("kraus"))
end

Py(ix::Instruction) = Py(ix.operator, ix.target)

function Py(ir::OpenQasmProgram)
    py_inputs = isnothing(ir.inputs) ? PythonCall.pybuiltins.None : pydict(ir.inputs)
    return braket[].ir.openqasm.program_v1.Program(source=pystr(ir.source), inputs=py_inputs)
end

function Py(ir::Program)
    bris = (isnothing(ir.basis_rotation_instructions) || isempty(ir.basis_rotation_instructions)) ? PythonCall.pybuiltins.None : pylist(ir.basis_rotation_instructions)
    res  = (isnothing(ir.results) || isempty(ir.results)) ? PythonCall.pybuiltins.None : pylist(ir.results)
    return braket[].ir.jaqcd.program_v1.Program(instructions=pylist(ir.instructions), results=res, basis_rotation_instructions=bris)
end

function Py(r::GateModelTaskResult)
    py_measurements  = (isnothing(r.measurements) || isempty(r.measurements)) ? PythonCall.pybuiltins.None : pylist(pylist(meas) for meas in r.measurements)
    py_probabilities = !isnothing(r.measurementProbabilities) ? pydict(Dict(pystr(k)=>v for (k,v) in r.measurementProbabilities)) : PythonCall.pybuiltins.None
    py_qubits        = !isnothing(r.measuredQubits) ? pylist(r.measuredQubits) : PythonCall.pybuiltins.None
    py_results       = pylist(Py(rtv) for rtv in r.resultTypes)
    py_task_mtd      = braket[].task_result.task_metadata_v1.TaskMetadata(id=pystr(r.taskMetadata.id), shots=Py(r.taskMetadata.shots), deviceId=pystr(r.taskMetadata.deviceId))
    # we create a placeholder OpenQASM3 Program here -- this is replaced by the "true" action in the Python wrapper
    dummy_action     = braket[].ir.openqasm.program_v1.Program(source=pystr(""))
    py_addl_mtd      = braket[].task_result.additional_metadata.AdditionalMetadata(action=dummy_action)
    return braket[].task_result.GateModelTaskResult(measurements=py_measurements, measurementProbabilities=py_probabilities, resultTypes=py_results, measuredQubits=py_qubits, taskMetadata=py_task_mtd, additionalMetadata=py_addl_mtd)
end
