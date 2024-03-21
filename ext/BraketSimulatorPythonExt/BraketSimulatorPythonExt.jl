module BraketSimulatorPythonExt

using PrecompileTools: @setup_workload, @compile_workload

using BraketSimulator, BraketSimulator.Braket, BraketSimulator.Braket.JSON3, PythonCall, BraketSimulator.Dates

import BraketSimulator.Braket:
    LocalSimulator,
    qubit_count,
    _run_internal,
    Instruction,
    Observables,
    AbstractProgramResult,
    ResultTypeValue,
    format_result,
    LocalQuantumTask,
    LocalQuantumTaskBatch,
    GateModelQuantumTaskResult,
    GateModelTaskResult,
    Program,
    Gate,
    AngledGate,
    AbstractIR,
    AbstractProgram,
    IRObservable
import BraketSimulator:
    AbstractSimulator,
    classical_shadow,
    AbstractStateVector,
    apply_gate!,
    get_amps_and_qubits,
    pad_bits,
    flip_bits,
    flip_bit,
    DoubleExcitation,
    SingleExcitation,
    MultiRZ

#const pennylane = Ref{Py}()
const numpy     = Ref{Py}()
const braket    = Ref{Py}()

include("translation.jl")

function __init__()
    # must set these when this code is actually loaded
    braket[]    = pyimport("braket")
    #pennylane[] = pyimport("pennylane")
    numpy[]     = pyimport("numpy")
    PythonCall.pyconvert_add_rule("braket.schema_common.schema_header:BraketSchemaHeader", Braket.braketSchemaHeader, jl_convert)
    PythonCall.pyconvert_add_rule("braket.circuits.circuit:Circuit", BraketSimulator.Braket.IR.Program, jl_convert_circuit)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.program_v1:Program", BraketSimulator.Braket.IR.Program, jl_convert_program)
    PythonCall.pyconvert_add_rule("braket.circuits.instruction:Instruction", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CNot", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Kraus", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDephasing", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDepolarizing", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:X", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift10", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift00", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:MultiQubitPauliChannel", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ti", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CV", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:StartVerbatimBox", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ECR", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CSwap", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ry", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CY", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CCNot", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PauliChannel", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:I", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Unitary", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Z", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Si", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift01", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:AmplitudeDamping", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PSwap", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:BitFlip", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseDamping", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rz", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:GeneralizedAmplitudeDamping", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseShift", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:V", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XX", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Y", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ZZ", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Swap", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ISwap", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:H", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseFlip", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:S", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Depolarizing", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rx", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:YY", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:EndVerbatimBox", EndVerbatimBox, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:T", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CZ", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XY", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Vi", Instruction, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CNot", CNot, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Kraus", Kraus, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDephasing", TwoQubitDephasing, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDepolarizing", TwoQubitDepolarizing, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:X", X, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift10", CPhaseShift10, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift00", CPhaseShift00, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:MultiQubitPauliChannel", MultiQubitPauliChannel, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ti", Ti, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CV", CV, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:StartVerbatimBox", StartVerbatimBox, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ECR", ECR, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CSwap", CSwap, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ry", Ry, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CY", CY, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CCNot", CCNot, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PauliChannel", PauliChannel, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:I", I, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Unitary", Unitary, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Z", Z, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Si", Si, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift01", CPhaseShift01, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:AmplitudeDamping", AmplitudeDamping, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PSwap", PSwap, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:BitFlip", BitFlip, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseDamping", PhaseDamping, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rz", Rz, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:GeneralizedAmplitudeDamping", GeneralizedAmplitudeDamping, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseShift", PhaseShift, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:V", V, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XX", XX, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Y", Y, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ZZ", ZZ, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Swap", Swap, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ISwap", ISwap, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:H", H, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift", CPhaseShift, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseFlip", PhaseFlip, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:S", S, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Depolarizing", Depolarizing, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rx", Rx, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:YY", YY, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:EndVerbatimBox", EndVerbatimBox, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:T", T, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CZ", CZ, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XY", XY, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Vi", Vi, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Sample", Sample, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Expectation", Expectation, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Probability", Probability, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:StateVector", Braket.IR.StateVector, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Amplitude", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Expectation", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Probability", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Sample", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:StateVector", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:DensityMatrix", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Variance", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:AdjointGradient", AbstractProgramResult, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:DensityMatrix", DensityMatrix, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Amplitude", Amplitude, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:AdjointGradient", AdjointGradient, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.shared_models:CompilerDirective", CompilerDirective, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.openqasm.program_v1:Program", OpenQasmProgram, jl_convert)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Variance", Variance, jl_convert)
    #=PythonCall.pyconvert_add_rule(
        "pennylane.ops.op_math:Adjoint",
        Instruction,
        pennylane_convert_Adjoint,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliX",
        Instruction,
        pennylane_convert_X,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliY",
        Instruction,
        pennylane_convert_Y,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliZ",
        Instruction,
        pennylane_convert_Z,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.identity:Identity",
        Instruction,
        pennylane_convert_I,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:Hadamard",
        Instruction,
        pennylane_convert_H,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:S",
        Instruction,
        pennylane_convert_S,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:SX",
        Instruction,
        pennylane_convert_V,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:T",
        Instruction,
        pennylane_convert_T,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:CNOT",
        Instruction,
        pennylane_convert_CNOT,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:CY",
        Instruction,
        pennylane_convert_CY,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:CZ",
        Instruction,
        pennylane_convert_CZ,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:Toffoli",
        Instruction,
        pennylane_convert_Toffoli,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:SWAP",
        Instruction,
        pennylane_convert_SWAP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:ISWAP",
        Instruction,
        pennylane_convert_ISWAP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:PSWAP",
        Instruction,
        pennylane_convert_PSWAP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:CSWAP",
        Instruction,
        pennylane_convert_CSWAP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:ECR",
        Instruction,
        pennylane_convert_ECR,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_single_qubit:RX",
        Instruction,
        pennylane_convert_RX,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_single_qubit:RY",
        Instruction,
        pennylane_convert_RY,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_single_qubit:RZ",
        Instruction,
        pennylane_convert_RZ,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_single_qubit:PhaseShift",
        Instruction,
        pennylane_convert_PhaseShift,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_controlled:ControlledPhaseShift",
        Instruction,
        pennylane_convert_ControlledPhaseShift,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_controlled:CPhaseShift00",
        Instruction,
        pennylane_convert_CPhaseShift00,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_controlled:CPhaseShift01",
        Instruction,
        pennylane_convert_CPhaseShift01,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_controlled:CPhaseShift10",
        Instruction,
        pennylane_convert_CPhaseShift10,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:IsingXX",
        Instruction,
        pennylane_convert_IsingXX,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:IsingXY",
        Instruction,
        pennylane_convert_IsingXY,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:IsingYY",
        Instruction,
        pennylane_convert_IsingYY,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:IsingZZ",
        Instruction,
        pennylane_convert_IsingZZ,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.parametric_ops_multi_qubit:MultiRZ",
        Instruction,
        pennylane_convert_MultiRZ,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.qchem_ops:DoubleExcitation",
        Instruction,
        pennylane_convert_DoubleExcitation,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.qchem_ops:SingleExcitation",
        Instruction,
        pennylane_convert_SingleExcitation,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.numpy.tensor:tensor",
        Union{Float64,FreeParameter},
        pennylane_convert_tensor,
    )
    PythonCall.pyconvert_add_rule(
        "builtins:list",
        Union{Float64,FreeParameter},
        pennylane_convert_parameters,
    )
    PythonCall.pyconvert_add_rule(
        "builtins:list",
        Union{Dict{String,Float64},Vector{Dict{String,Float64}}},
        pennylane_convert_inputs,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliX",
        Observables.Observable,
        pennylane_convert_X,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliY",
        Observables.Observable,
        pennylane_convert_Y,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliZ",
        Observables.Observable,
        pennylane_convert_Z,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.identity:Identity",
        Observables.Observable,
        pennylane_convert_I,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:Hadamard",
        Observables.Observable,
        pennylane_convert_H,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.observables:Hermitian",
        Observables.Observable,
        pennylane_convert_Hermitian,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.operation:Tensor",
        Observables.Observable,
        pennylane_convert_Tensor,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliX",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_X,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliY",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_Y,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:PauliZ",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_Z,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.identity:Identity",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_I,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.non_parametric_ops:Hadamard",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_H,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.ops.qubit.observables:Hermitian",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_Hermitian,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.operation:Tensor",
        Tuple{IRObservable,Vector{Int}},
        pennylane_convert_Tensor,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.measurements.expval:ExpectationMP",
        BraketSimulator.Braket.IR.AbstractProgramResult,
        pennylane_convert_ExpectationMP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.measurements.var:VarianceMP",
        BraketSimulator.Braket.IR.AbstractProgramResult,
        pennylane_convert_VarianceMP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.measurements.sample:SampleMP",
        BraketSimulator.Braket.IR.AbstractProgramResult,
        pennylane_convert_SampleMP,
    )
    PythonCall.pyconvert_add_rule(
        "pennylane.tape.qscript:QuantumScript",
        BraketSimulator.Braket.IR.Program,
        pennylane_convert_QuantumScript,
    )=#
end
#BraketSimulator.Braket.qubit_count(o::Py) =
#    pyisinstance(o, pennylane.tape.QuantumTape) ? length(o.wires) : o.qubit_count

function classical_shadow(d::LocalSimulator, obs_qubits, circuit, shots::Int, seed::Int)
    raw_jl_spec = _translate_from_python(circuit, d._delegate)
    PythonCall.GC.disable()
    jl_spec = pyconvert(Program, circuit)
    shadow = classical_shadow(d, pyconvert(Vector{Int}, obs_qubits), jl_spec, shots, seed)
    PythonCall.GC.enable()
    return shadow
end

function classical_shadow(
    d::LocalSimulator,
    obs_qubits,
    circuits::PyList{Any},
    shots::Int,
    seed::PyList,
)
    jl_obs_qubits = pyconvert(Vector{Vector{Int}}, obs_qubits)
    jl_seed = pyconvert(Vector{Int}, seed)
    jl_specs = [pyconvert(Program, circuit) for circuit in circuits]
    PythonCall.GC.disable()
    shadow = classical_shadow(d, jl_obs_qubits, jl_specs, shots, jl_seed)
    PythonCall.GC.enable()
    return shadow
end

function (d::AbstractSimulator)(
    task_specs::Union{PyList{Any},NTuple{N,PyIterable}, Py},
    args...;
    input::Union{PyList{Any},PyDict{Any,Any},Py}=PyDict{Any,Any}(),
    kwargs...,
) where {N}
    # handle inputs
    # fix me to use OpenQASM if needed
    jl_specs  = [] 
    jl_inputs = nothing
    shots = args[end]
    stats = @timed begin
        if input isa PyDict{Any,Any}
            jl_inputs = pyconvert(Dict{String,Float64}, input)
        else
            jl_inputs = [pyconvert(Dict{String,Float64}, py_inputs) for py_inputs in input]
        end
        s_ix = 1
        for spec in task_specs
            if pyisinstance(spec, braket[].ir.openqasm.Program)
                push!(jl_specs, pyconvert(OpenQasmProgram, spec))
            else
                push!(jl_specs, pyconvert(Program, spec))
            end
            s_ix += 1
        end
        task_specs = nothing
        input = nothing
    end
    @debug "Time for conversion of specs and inputs: $(stats.time)."
    PythonCall.GC.disable()
    if length(jl_specs) == 1
        r = d(jl_specs[1], args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
    else
        r = d(jl_specs, args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
    end
    PythonCall.GC.enable()
    return r
end

function Py(r::Braket.IR.Sample)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    raw_obs = map(r.observable) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    py_obs = pylist(raw_obs)
    return braket[].ir.jaqcd.results.Sample(targets=py_targets, observable=py_obs, type=pystr("sample"))
end

function Py(r::Braket.IR.Expectation)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    raw_obs = map(r.observable) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    py_obs = pylist(raw_obs)
    return braket[].ir.jaqcd.results.Expectation(targets=py_targets, observable=py_obs, type=pystr("expectation"))
end

function Py(r::Braket.IR.Variance)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    raw_obs = map(r.observable) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    py_obs = pylist(raw_obs)
    return braket[].ir.jaqcd.results.Variance(targets=py_targets, observable=py_obs, type=pystr("variance"))
end

function Py(r::Braket.IR.Amplitude)
    return braket[].ir.jaqcd.results.Amplitude(states=pylist(pystr(s) for s in r.states), type=pystr("amplitude"))
end

function Py(r::Braket.IR.StateVector)
    return braket[].ir.jaqcd.results.StateVector(type=pystr("statevector"))
end

function Py(r::Braket.IR.DensityMatrix)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.DensityMatrix(targets=py_targets, type=pystr("densitymatrix"))
end

function Py(r::Braket.IR.Probability)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.Probability(targets=py_targets, type=pystr("probability"))
end

function Py(r::Braket.IR.AdjointGradient)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.AdjointGradient(targets=py_targets, observable=pylist(r.observable), parameters=pylist(pystr(p) for p in r.parameters), type=pystr("adjoint_gradient"))
end

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
Py(op::Braket.XX, targets) = braket[].ir.jaqcd.instructions.XX(target=pylist(targets), angle=op.angle[1], type=pystr("xx"))
Py(op::Braket.YY, targets) = braket[].ir.jaqcd.instructions.YY(target=pylist(targets), angle=op.angle[1], type=pystr("yy"))
Py(op::Braket.ZZ, targets) = braket[].ir.jaqcd.instructions.ZZ(target=pylist(targets), angle=op.angle[1], type=pystr("zz"))
Py(op::Braket.XY, targets) = braket[].ir.jaqcd.instructions.XY(target=pylist(targets), angle=op.angle[1], type=pystr("xy"))
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
Py(op::Braket.CPhaseShift01, targets) = braket[].ir.jaqcd.instructions.CPhaseShift00(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift01"))
Py(op::Braket.CPhaseShift10, targets) = braket[].ir.jaqcd.instructions.CPhaseShift00(control=targets[1], target=targets[2], angle=op.angle[1], type=pystr("cphaseshift10"))
Py(op::Braket.ECR, targets) = braket[].ir.jaqcd.instructions.ECR(targets=pylist(targets), type=pystr("ecr"))
Py(op::Braket.BitFlip, targets) = braket[].ir.jaqcd.instructions.BitFlip(target=targets[1], probability=op.probability, type=pystr("bit_flip"))
Py(op::Braket.PhaseFlip, targets) = braket[].ir.jaqcd.instructions.PhaseFlip(target=targets[1], probability=op.probability, type=pystr("phase_flip"))
Py(op::Braket.PauliChannel, targets) = braket[].ir.jaqcd.instructions.PauliChannel(target=targets[1], probX=op.probX, probY=op.probY, probZ=op.probZ, type=pystr("pauli_channel"))
Py(op::Braket.MultiQubitPauliChannel, targets) = braket[].ir.jaqcd.instructions.MultiQubitPauliChannel(target=pylist(targets), probabilities=pydict(op.probabilities), type=pystr("multi_qubit_pauli_channel"))
Py(op::Braket.Depolarizing, targets) = braket[].ir.jaqcd.instructions.Depolarizing(target=targets[1], probability=op.probability, type=pystr("depolarizing"))
Py(op::Braket.TwoQubitDepolarizing, targets) = braket[].ir.jaqcd.instructions.TwoQubitDepolarizing(target1=targets[1], target2=targets[2], probability=op.probability, type=pystr("two_qubit_depolarizing"))
Py(op::Braket.TwoQubitDephasing, targets) = braket[].ir.jaqcd.instructions.TwoQubitDephasing(target1=targets[1], target2=targets[2], probability=op.probability, type=pystr("two_qubit_dephasing"))
Py(op::Braket.AmplitudeDamping, targets) = braket[].ir.jaqcd.instructions.AmplitudeDamping(target=targets[1], gamma=op.gamma, type=pystr("amplitude_damping"))
Py(op::Braket.PhaseDamping, targets) = braket[].ir.jaqcd.instructions.PhaseDamping(target=targets[1], gamma=op.gamma, type=pystr("phase_damping"))
Py(op::Braket.GeneralizedAmplitudeDamping, targets) = braket[].ir.jaqcd.instructions.GeneralizedAmplitudeDamping(target=targets[1], gamma=op.gamma, probability=op.probability, type=pystr("generalized_amplitude_damping"))
Py(op::Braket.StartVerbatimBox) = braket[].ir.jaqcd.instructions.StartVerbatimBox(type=pystr("start_verbatim_box"))
Py(op::Braket.EndVerbatimBox) = braket[].ir.jaqcd.instructions.EndVerbatimBox(type=pystr("end_verbatim_box"))
function Py(op::Braket.Unitary, targets)
    raw_mat = Braket.complex_matrix_to_ir(op.matrix)
    py_mat = pylist(pylist(pylist(v__) for v__ in v_) for v_ in raw_mat)
    return braket[].ir.jaqcd.instructions.Unitary(targets=pylist(targets), matrix=py_mat, type=pystr("unitary"))
end
function Py(op::Braket.Kraus, targets)
    raw_mats = map(Braket.complex_matrix_to_ir, op.matrices)
    py_mat = pylist(pylist(pylist(pylist(v__) for v__ in v_) for v_ in raw_mat) for raw_mat in raw_mats)
    return braket[].ir.jaqcd.instructions.Kraus(targets=pylist(targets), matrices=py_mats, type=pystr("kraus"))
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
    py_measurements  = isempty(r.measurements) ? PythonCall.pybuiltins.None : pylist(pylist(meas) for meas in r.measurements)
    py_probabilities = !isnothing(r.measurementProbabilities) ? pydict(Dict(pystr(k)=>v for (k,v) in r.measurementProbabilities)) : PythonCall.pybuiltins.None
    py_qubits = !isnothing(r.measuredQubits) ? pylist(r.measuredQubits) : PythonCall.pybuiltins.None
    py_results = pylist(r.resultTypes)
    py_task_mtd = braket[].task_result.task_metadata_v1.TaskMetadata(id=pystr(r.taskMetadata.id), shots=Py(r.taskMetadata.shots), deviceId=pystr(r.taskMetadata.deviceId))
    py_addl_mtd = braket[].task_result.additional_metadata.AdditionalMetadata(action=Py(r.additionalMetadata.action))
    return braket[].task_result.GateModelTaskResult(measurements=py_measurements, measurementProbabilities=py_probabilities, resultTypes=py_results, measuredQubits=py_qubits, taskMetadata=py_task_mtd, additionalMetadata=py_addl_mtd)
end
# PL specific -- some way we can dispatch here?
function Py(r::GateModelQuantumTaskResult)
    return pylist([numpy[].array(v).squeeze() for v in r.values])
end

@setup_workload begin
    @compile_workload begin
        pybraket = pyimport("braket.circuits")
        c = pybraket.Circuit()
        c.h(0)
        n_qubits = 10
        for q in 1:n_qubits-1
            c.cnot(0, q)
        end
        c.state_vector()
        svs = StateVectorSimulator(n_qubits, 0)
        svs([c], shots=0)
    end
end

end
