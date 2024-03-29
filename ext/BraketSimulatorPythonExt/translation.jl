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
jl_convert_kraus(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Kraus([convert_ir_matrix(m) for m in x.matrices]), [t for t in x.targets]))

jl_convert_unitary(::Type{Unitary}, x::Py)::Unitary = PythonCall.pyconvert_return(Unitary(convert_ir_matrix(x.matrix)))
jl_convert_unitary(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(Unitary(convert_ir_matrix(x.matrix)), [t for t in x.targets]))

jl_convert_cswap(::Type{CSwap}, x::Py)::CSwap = PythonCall.pyconvert_return(CSwap())
jl_convert_cswap(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(CSwap(), [t for t in x.targets]))
jl_convert_ccnot(::Type{CCNot}, x::Py)::CCNot = PythonCall.pyconvert_return(CCNot())
jl_convert_ccnot(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(CCNot(), vcat([c for c in x.controls], x.target)))

jl_convert_paulichannel(::Type{PauliChannel}, x::Py)::PauliChannel = PythonCall.pyconvert_return(PauliChannel(x.probX, x.probY, x.probZ))
jl_convert_paulichannel(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(PauliChannel(x.probX, x.probY, x.probZ), x.target))

jl_convert_generalizedampdamp(::Type{GeneralizedAmplitudeDamping}, x::Py)::GeneralizedAmplitudeDamping = PythonCall.pyconvert_return(GeneralizedAmplitudeDamping(x.probability, x.gamma))
jl_convert_generalizedampdamp(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction(GeneralizedAmplitudeDamping(x.probability, x.gamma), x.target))

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
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(x.probability))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.probability)), pyconvert(Int, x.target)))
    end
end
for (g, fn) in ((:AmplitudeDamping, :jl_convert_amplitudedamp),
                (:PhaseDamping, :jl_convert_phasedamp),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(x.gamma))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.gamma)), pyconvert(Int, x.target)))
    end
end
for (g, fn) in ((:TwoQubitDepolarizing, :jl_convert_twoqubitdepo),
                (:TwoQubitDephasing, :jl_convert_twoqubitdeph),
               )
    @eval begin
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(x.probability))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.probability)), [pyconvert(Int, t) for t in x.targets]))
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
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.angle)), pyconvert(Int, x.target)))
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
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.angle)), [pyconvert(Int, t) for t in x.targets]))
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
        $fn(::Type{$g}, x::Py)::$g = PythonCall.pyconvert_return($g(pyconvert(Float64, x.angle)))
        $fn(::Type{Instruction}, x::Py)::Instruction = PythonCall.pyconvert_return(Instruction($g(pyconvert(Float64, x.angle)), [pyconvert(Int, x.control), pyconvert(Int, x.target)]))
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
#=pennylane_convert_tensor(::Type{Float64}, x::Py) =
    PythonCall.pyconvert_return(pyconvert(Float64, x.numpy()))
pennylane_convert_parameters(::Type{Float64}, x::Py) =
    PythonCall.pyconvert_return(pyconvert(Float64, x[0]))
pennylane_convert_parameters(::Type{FreeParameter}, x::Py) =
    PythonCall.pyconvert_return(FreeParameter(pyconvert(String, x[0])))

pennylane_convert_inputs(::Type{Dict{String,Float64}}, x::Py) =
    PythonCall.pyconvert_return(pyconvert(Dict{String,Float64}, x))
pennylane_convert_inputs(::Type{Vector{Dict{String,Float64}}}, x::Py) =
    PythonCall.pyconvert_return([pyconvert(Dict{String,Float64}, x_) for x_ in x])
for (conv_fn, jl_typ) in (
    (:pennylane_convert_X, :X),
    (:pennylane_convert_Y, :Y),
    (:pennylane_convert_Z, :Z),
    (:pennylane_convert_I, :I),
    (:pennylane_convert_H, :H),
    (:pennylane_convert_S, :S),
    (:pennylane_convert_T, :T),
    (:pennylane_convert_V, :V),
)
    @eval begin
        function $conv_fn(::Type{Instruction}, x::Py)
            return PythonCall.pyconvert_return(
                Instruction($jl_typ(), pyconvert(Int, x.wires[0])),
            )
        end
    end
end

for (conv_fn, jl_typ) in (
    (:pennylane_convert_RX, :Rx),
    (:pennylane_convert_RY, :Ry),
    (:pennylane_convert_RZ, :Rz),
    (:pennylane_convert_PhaseShift, :PhaseShift),
)
    @eval begin
        function $conv_fn(::Type{Instruction}, x::Py)
            angle = pyconvert(Union{Float64,FreeParameter}, x.parameters)
            return PythonCall.pyconvert_return(
                Instruction($jl_typ(angle), pyconvert(Int, x.wires[0])),
            )
        end
    end
end

for (conv_fn, jl_typ) in (
    (:pennylane_convert_IsingXX, :XX),
    (:pennylane_convert_IsingYY, :YY),
    (:pennylane_convert_IsingZZ, :ZZ),
    (:pennylane_convert_IsingXY, :XY),
    (:pennylane_convert_PSWAP, :PSwap),
    (:pennylane_convert_ControlledPhaseShift, :CPhaseShift),
    (:pennylane_convert_CPhaseShift00, :CPhaseShift00),
    (:pennylane_convert_CPhaseShift01, :CPhaseShift01),
    (:pennylane_convert_CPhaseShift10, :CPhaseShift10),
    (:pennylane_convert_DoubleExcitation, :DoubleExcitation),
    (:pennylane_convert_SingleExcitation, :SingleExcitation),
)
    @eval begin
        function $conv_fn(::Type{Instruction}, x::Py)
            angle = pyconvert(Union{Float64,FreeParameter}, x.parameters)
            return PythonCall.pyconvert_return(
                Instruction($jl_typ(angle), pyconvert(Vector{Int}, x.wires)),
            )
        end
    end
end
function pennylane_convert_MultiRZ(::Type{Instruction}, x::Py)
    angle = pyconvert(Union{Float64,FreeParameter}, x.parameters)
    return PythonCall.pyconvert_return(
        Instruction(MultiRZ((angle,),), pyconvert(Vector{Int}, x.wires)),
    )
end

for (conv_fn, jl_typ) in (
    (:pennylane_convert_CNOT, :CNot),
    (:pennylane_convert_CY, :CY),
    (:pennylane_convert_CZ, :CZ),
    (:pennylane_convert_SWAP, :Swap),
    (:pennylane_convert_ISWAP, :ISwap),
    (:pennylane_convert_CSWAP, :CSwap),
    (:pennylane_convert_ECR, :ECR),
    (:pennylane_convert_Toffoli, :CCNot),
)
    @eval begin
        function $conv_fn(::Type{Instruction}, x::Py)
            return PythonCall.pyconvert_return(
                Instruction($jl_typ(), pyconvert(Vector{Int}, x.wires)),
            )
        end
    end
end
function pennylane_convert_QubitUnitary(::Type{Instruction}, x::Py)
    mat = pyconvert(Matrix{ComplexF64}, x.parameters[0])
    return PythonCall.pyconvert_return(Instruction(Unitary(mat), pyconvert(Int, x.wires)))
end

for (typ, adj_typ) in ((:S, :Si), (:T, :Ti), (:V, :Vi))
    @eval begin
        adjoint_type(::Type{$typ}) = $adj_typ
        adjoint_type(g::$typ) = $adj_typ()
    end
end

function pennylane_convert_Adjoint(::Type{Instruction}, x::Py)
    un_adjointed_instruction = pyconvert(Instruction, x.base)
    raw_gate = un_adjointed_instruction.operator
    return PythonCall.pyconvert_return(
        Instruction(adjoint_type(raw_gate), un_adjointed_instruction.target),
    )
end

function _translate_parameters(py_params, parameter_names::Vector{String}, ::Val{true})
    isempty(py_params) && return Float64[]
    param_names = isempty(parameter_names) ? fill("", length(py_params)) : parameter_names
    length(param_names) != length(py_params) && throw(
        ErrorException(
            "Parameter names list must be equal to number of operation parameters",
        ),
    )
    parameters = map(zip(param_names, py_params)) do (param_name, param)
        # PennyLane passes any non-keyword argument in the operation.parameters list.
        # In some cases, like the unitary gate or qml.QubitChannel (Kraus noise), these
        # parameter can be matrices. Braket only supports parameterization of numeric parameters
        # (so far, these are all angle parameters), so non-numeric parameters are handled
        # separately.
        param_name != "" && return BraketStateVector.Braket.FreeParameter(param_name)
        pyisinstance(param, pennylane.numpy.tensor) &&
            return pyconvert(Array, param.numpy())
        return pyconvert(Float64, param)
    end
    return parameters
end

for (conv_fn, jl_typ, str) in (
    (:pennylane_convert_X, :(Observables.X), "x"),
    (:pennylane_convert_Y, :(Observables.Y), "y"),
    (:pennylane_convert_Z, :(Observables.Z), "z"),
    (:pennylane_convert_I, :(Observables.I), "i"),
    (:pennylane_convert_H, :(Observables.H), "h"),
)
    @eval begin
        $conv_fn(::Type{Observables.Observable}, x::Py) =
            PythonCall.pyconvert_return($jl_typ())
        $conv_fn(::Type{Tuple{IRObservable,Vector{Int}}}, x::Py) =
            PythonCall.pyconvert_return(($str, pyconvert(Vector{Int}, x.wires)))
    end
end
function pennylane_convert_Hermitian(::Type{Observables.Observable}, o::Py)
    return PythonCall.pyconvert_return(
        BraketStateVector.Braket.Observables.HermitianObservable(
            pyconvert(Matrix{ComplexF64}, o.parameters[0]),
        ),
    )
end
function pennylane_convert_Hermitian(::Type{Tuple{IRObservable,Vector{Int}}}, o::Py)
    mat = BraketStateVector.Braket.complex_matrix_to_ir(
        pyconvert(Matrix{ComplexF64}, o.parameters[0]),
    )
    return PythonCall.pyconvert_return((mat, pyconvert(Vector{Int}, o.wires)))
end

function pennylane_convert_Tensor(::Type{Observables.Observable}, o::Py)
    return PythonCall.pyconvert_return(
        Observables.TensorProduct([pyconvert(Observables.Observable, o.obs)]),
    )
end
function pennylane_convert_Tensor(::Type{Tuple{IRObservable,Vector{Int}}}, o::Py)
    raw_obs = [pyconvert(Tuple{IRObservable,Vector{Int}}, f) for f in o.obs]
    tensor_ops = convert(IRObservable, reduce(vcat, [o[1] for o in raw_obs]))
    tensor_qubits = reduce(vcat, [o[2] for o in raw_obs])
    return PythonCall.pyconvert_return((tensor_ops, tensor_qubits))
end

for (ir_typ, conv_fn, braket_name) in (
    (
        :(BraketStateVector.Braket.IR.Expectation),
        :pennylane_convert_ExpectationMP,
        "expectation",
    ),
    (:(BraketStateVector.Braket.IR.Variance), :pennylane_convert_VarianceMP, "variance"),
    (:(BraketStateVector.Braket.IR.Sample), :pennylane_convert_SampleMP, "sample"),
)
    @eval begin
        function $conv_fn(::Type{AbstractProgramResult}, o::Py)
            ir_obs, ir_qubits = pyconvert(Tuple{IRObservable,Vector{Int}}, o.obs)
            return PythonCall.pyconvert_return($ir_typ(ir_obs, ir_qubits, $braket_name))
        end
    end
end

function pennylane_convert_QuantumScript(::Type{Program}, o)
    instructions = [pyconvert(Instruction, i) for i in o.operations]
    results_list = [pyconvert(AbstractProgramResult, i) for i in o.measurements]
    instr_qubits = mapreduce(ix -> ix.target, union, instructions)
    result_qubits = mapreduce(
        ix -> hasproperty(ix, :targets) ? ix.targets : Set{Int}(),
        union,
        results_list,
        init = Set{Int}(),
    )
    all_qubits = union(result_qubits, instr_qubits)
    missing_qubits = union(
        setdiff(result_qubits, instr_qubits),
        setdiff(0:maximum(all_qubits), instr_qubits),
    )
    for q in missing_qubits
        push!(instructions, Instruction(Braket.I(), q))
    end
    prog = Program(
        BraketStateVector.Braket.header_dict[Program],
        instructions,
        results_list,
        [],
    )
    return PythonCall.pyconvert_return(prog)
end

function _translate_parameter_names(
    n_params::Int,
    param_index::Int,
    trainable_indices::Set{Int},
    use_unique_parameters::Bool,
    ::Val{false},
)
    n_params == 0 && return String[], param_index
    parameter_names = fill("", n_params)
    ix = 1
    for p = 1:n_params
        if param_index ∈ trainable_indices || use_unique_parameters
            parameter_names[ix] = "p_$param_index"
            ix += 1
        end
        param_index += 1
    end
    return parameter_names, param_index
end

function _translate_parameter_names(
    n_params::Int,
    param_index::Int,
    trainable_indices::Set{Int},
    use_unique_parameters::Bool,
    ::Val{true},
)
    return fill("", n_params), param_index + n_params
end=#
