module BraketSimulatorPythonExt

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

const pennylane = Ref{Py}()
const numpy     = Ref{Py}()
const braket    = Ref{Py}()

include("translation.jl")

function __init__()
    # must set these when this code is actually loaded
    braket[]    = pyimport("braket")
    pennylane[] = pyimport("pennylane")
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
    PythonCall.pyconvert_add_rule(
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
    )
end
BraketSimulator.Braket.qubit_count(o::Py) =
    pyisinstance(o, pennylane.tape.QuantumTape) ? length(o.wires) : o.qubit_count

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
    return braket[].ir.jaqcd.results.Sample(targets=pylist(r.targets), observable=pylist(r.observable))
end

function Py(r::Braket.IR.Expectation)
    return braket[].ir.jaqcd.results.Expectation(targets=pylist(r.targets), observable=pylist(r.observable))
end

function Py(r::Braket.IR.Variance)
    return braket[].ir.jaqcd.results.Variance(targets=pylist(r.targets), observable=pylist(r.observable))
end

function Py(r::Braket.IR.Amplitude)
    return braket[].ir.jaqcd.results.Amplitude(states=pylist(pystr(s) for s in r.states))
end

function Py(r::Braket.IR.StateVector)
    return braket[].ir.jaqcd.results.StateVector()
end

function Py(r::Braket.IR.DensityMatrix)
    return braket[].ir.jaqcd.results.DensityMatrix(targets=pylist(r.targets))
end

function Py(r::Braket.IR.Probability)
    return braket[].ir.jaqcd.results.Probability(targets=pylist(r.targets))
end

function Py(r::Braket.IR.AdjointGradient)
    return braket[].ir.jaqcd.results.AdjointGradient(targets=pylist(r.targets), observable=pylist(r.observable), parameters=pylist(pystr(p) for p in r.parameters))
end

function Py(rt::Braket.ResultTypeValue)
    py_typ = Py(rt.type) 
    py_val = if rt.value isa Dict
        pydict(rt.value)
    elseif rt.value isa Float64
        rt.value
    else
        pylist(rt.value)
    end
    return braket[].task_result.gate_model_task_result_v1.ResultTypeValue(type=py_typ, value=py_val)
end

function Py(ir::OpenQasmProgram)
    py_inputs = isnothing(ir.inputs) ? PythonCall.pybuiltins.None : pydict(ir.inputs)
    return braket[].ir.openqasm.program_v1.Program(source=pystr(ir.source), inputs=py_inputs)
end

function Py(r::GateModelTaskResult)
    py_measurements  = pylist( pylist(Py(m) for m in meas) for meas in r.measurements)
    py_probabilities = !isnothing(r.measurementProbabilities) ? pydict(r.measurementProbabilities) : PythonCall.pybuiltins.None
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

end
