module BraketSimulatorPythonExt

using PrecompileTools: @setup_workload, @compile_workload, @recompile_invalidations

using BraketSimulator, BraketSimulator.Braket, PythonCall, BraketSimulator.Dates

import BraketSimulator.Braket:
    LocalSimulator,
    qubit_count,
    Instruction,
    Observables,
    AbstractProgramResult,
    ResultTypeValue,
    format_result,
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
    simulate,
    parse_program,
    DoubleExcitation,
    SingleExcitation,
    Control,
    MultiQubitPhaseShift,
    MultiRZ

const numpy     = Ref{Py}()
const braket    = Ref{Py}()
const sympy     = Ref{Py}()

include("translation.jl")

function __init__()
    # must set these when this code is actually loaded
    braket[]    = pyimport("braket")
    numpy[]     = pyimport("numpy")
    sympy[]     = pyimport("sympy")
    PythonCall.pyconvert_add_rule("sympy.core.mul:Mul", Union{Float64, Braket.FreeParameter}, jl_convert_sympy_Mul)
    PythonCall.pyconvert_add_rule("sympy.core.numbers:Pi", Float64, jl_convert_sympy_Pi)
    PythonCall.pyconvert_add_rule("sympy.core.numbers:Exp1", Float64, jl_convert_sympy_E)
    PythonCall.pyconvert_add_rule("braket.schema_common.schema_header:BraketSchemaHeader", Braket.braketSchemaHeader, jl_convert_bsh)
    PythonCall.pyconvert_add_rule("braket.circuits.circuit:Circuit", BraketSimulator.Braket.IR.Program, jl_convert_circuit)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.program_v1:Program", BraketSimulator.Braket.IR.Program, jl_convert_program)
    PythonCall.pyconvert_add_rule("braket.circuits.instruction:Instruction", Instruction, jl_convert_ix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.openqasm.circuit:Circuit", Circuit, jl_convert_sim_circuit)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Identity", Instruction, jl_convert_sim_identity)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Hadamard", Instruction, jl_convert_sim_hadamard)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:PauliX", Instruction, jl_convert_sim_paulix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:PauliY", Instruction, jl_convert_sim_pauliy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:PauliZ", Instruction, jl_convert_sim_pauliz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CV", Instruction, jl_convert_sim_cv)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CX", Instruction, jl_convert_sim_cx)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CY", Instruction, jl_convert_sim_cy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CZ", Instruction, jl_convert_sim_cz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:ECR", Instruction, jl_convert_sim_ecr)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:MS", Instruction, jl_convert_sim_ms)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:GPi", Instruction, jl_convert_sim_gpi)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:GPi2", Instruction, jl_convert_sim_gpi2)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:S", Instruction, jl_convert_sim_s)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Si", Instruction, jl_convert_sim_si)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:V", Instruction, jl_convert_sim_v)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Vi", Instruction, jl_convert_sim_vi)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:T", Instruction, jl_convert_sim_t)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Ti", Instruction, jl_convert_sim_ti)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:PhaseShift", Instruction, jl_convert_sim_phaseshift)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CPhaseShift", Instruction, jl_convert_sim_cphaseshift)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CPhaseShift00", Instruction, jl_convert_sim_cphaseshift00)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CPhaseShift01", Instruction, jl_convert_sim_cphaseshift01)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CPhaseShift10", Instruction, jl_convert_sim_cphaseshift10)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:RotX", Instruction, jl_convert_sim_rx)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:RotY", Instruction, jl_convert_sim_ry)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:RotZ", Instruction, jl_convert_sim_rz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Swap", Instruction, jl_convert_sim_swap)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:ISwap", Instruction, jl_convert_sim_iswap)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:PSwap", Instruction, jl_convert_sim_pswap)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CSwap", Instruction, jl_convert_sim_cswap)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:XY", Instruction, jl_convert_sim_xy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:XX", Instruction, jl_convert_sim_xx)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:YY", Instruction, jl_convert_sim_yy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:ZZ", Instruction, jl_convert_sim_zz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:CCNot", Instruction, jl_convert_sim_ccnot)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:Unitary", Instruction, jl_convert_sim_unitary)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:U", Instruction, jl_convert_sim_u)
    PythonCall.pyconvert_add_rule("braket.default_simulator.gate_operations:GPhase", Instruction, jl_convert_sim_gphase)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:BitFlip", Instruction, jl_convert_sim_bitflip)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:PhaseFlip", Instruction, jl_convert_sim_phaseflip)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:PauliChannel", Instruction, jl_convert_sim_paulichannel)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:Depolarizing", Instruction, jl_convert_sim_depolarizing)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:TwoQubitDepolarizing", Instruction, jl_convert_sim_twoqubitdepolarizing)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:TwoQubitDephasing", Instruction, jl_convert_sim_twoqubitdephasing)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:AmplitudeDamping", Instruction, jl_convert_sim_amplitudedamping)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:GeneralizedAmplitudeDamping", Instruction, jl_convert_sim_generalizedamplitudedamping)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:PhaseDamping", Instruction, jl_convert_sim_phasedamping)
    PythonCall.pyconvert_add_rule("braket.default_simulator.noise_operations:Kraus", Instruction, jl_convert_sim_kraus)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Expectation", Braket.Expectation, jl_convert_sim_expectation)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Variance", Braket.Variance, jl_convert_sim_variance)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Probability", Braket.Probability, jl_convert_sim_probability)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:StateVector", Braket.StateVector, jl_convert_sim_statevector)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:DensityMatrix", Braket.DensityMatrix, jl_convert_sim_densitymatrix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Amplitude", Braket.Amplitude, jl_convert_sim_amplitude)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Amplitude", Braket.Result, jl_convert_sim_amplitude)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Expectation", Braket.Result, jl_convert_sim_expectation)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Probability", Braket.Result, jl_convert_sim_probability)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:StateVector", Braket.Result, jl_convert_sim_statevector)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:DensityMatrix", Braket.Result, jl_convert_sim_densitymatrix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:Variance", Braket.Result, jl_convert_sim_variance)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:AdjointGradient", Braket.Result, jl_convert_sim_adjointgradient)
    PythonCall.pyconvert_add_rule("braket.default_simulator.result_types:AdjointGradient", Braket.AdjointGradient, jl_convert_sim_adjointgradient)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Identity", Braket.Observables.I, jl_convert_sim_identity)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Identity", Braket.Observables.Observable, jl_convert_sim_identity)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Hadamard", Braket.Observables.H, jl_convert_sim_hadamard)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Hadamard", Braket.Observables.Observable, jl_convert_sim_hadamard)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliX", Braket.Observables.X, jl_convert_sim_paulix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliX", Braket.Observables.Observable, jl_convert_sim_paulix)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliY", Braket.Observables.Y, jl_convert_sim_pauliy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliY", Braket.Observables.Observable, jl_convert_sim_pauliy)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliZ", Braket.Observables.Z, jl_convert_sim_pauliz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:PauliZ", Braket.Observables.Observable, jl_convert_sim_pauliz)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Hermitian", Braket.Observables.HermitianObservable, jl_convert_sim_hermitian)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:Hermitian", Braket.Observables.Observable, jl_convert_sim_hermitian)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:TensorProduct", Braket.Observables.TensorProduct, jl_convert_sim_tensorproduct)
    PythonCall.pyconvert_add_rule("braket.default_simulator.observables:TensorProduct", Braket.Observables.Observable, jl_convert_sim_tensorproduct)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CNot", Instruction, jl_convert_cnot)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Kraus", Instruction, jl_convert_kraus)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDephasing", Instruction, jl_convert_twoqubitdeph)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDepolarizing", Instruction, jl_convert_twoqubitdepo)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:X", Instruction, jl_convert_x)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift10", Instruction, jl_convert_cphaseshift10)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift00", Instruction, jl_convert_cphaseshift00)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:MultiQubitPauliChannel", Instruction, jl_convert_multi_qubit_pauli_channel)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ti", Instruction, jl_convert_ti)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CV", Instruction, jl_convert_cv)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:StartVerbatimBox", Instruction, jl_convert_startverbatim)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ECR", Instruction, jl_convert_ecr)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CSwap", Instruction, jl_convert_cswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ry", Instruction, jl_convert_ry)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CY", Instruction, jl_convert_cy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CCNot", Instruction, jl_convert_ccnot)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PauliChannel", Instruction, jl_convert_paulichannel)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:I", Instruction, jl_convert_i)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Unitary", Instruction, jl_convert_unitary)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Z", Instruction, jl_convert_z)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Si", Instruction, jl_convert_si)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift01", Instruction, jl_convert_cphaseshift01)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:AmplitudeDamping", Instruction, jl_convert_amplitudedamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PSwap", Instruction, jl_convert_pswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:BitFlip", Instruction, jl_convert_bitflip)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseDamping", Instruction, jl_convert_phasedamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rz", Instruction, jl_convert_rz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:GeneralizedAmplitudeDamping", Instruction, jl_convert_generalizedampdamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseShift", Instruction, jl_convert_phaseshift)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:V", Instruction, jl_convert_v)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XX", Instruction, jl_convert_xx)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Y", Instruction, jl_convert_y)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ZZ", Instruction, jl_convert_zz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Swap", Instruction, jl_convert_swap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ISwap", Instruction, jl_convert_iswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:H", Instruction, jl_convert_h)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift", Instruction, jl_convert_cphaseshift)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseFlip", Instruction, jl_convert_phaseflip)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:S", Instruction, jl_convert_s)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Depolarizing", Instruction, jl_convert_depo)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rx", Instruction, jl_convert_rx)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:YY", Instruction, jl_convert_yy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:EndVerbatimBox", Instruction, jl_convert_endverbatim)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:T", Instruction, jl_convert_t)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CZ", Instruction, jl_convert_cz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XY", Instruction, jl_convert_xy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Vi", Instruction, jl_convert_vi)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CNot", CNot, jl_convert_cnot)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Kraus", Kraus, jl_convert_kraus)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDephasing", TwoQubitDephasing, jl_convert_twoqubitdeph)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:TwoQubitDepolarizing", TwoQubitDepolarizing, jl_convert_twoqubitdepo)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:X", X, jl_convert_x)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift10", CPhaseShift10, jl_convert_cphaseshift10)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift00", CPhaseShift00, jl_convert_cphaseshift00)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:MultiQubitPauliChannel", MultiQubitPauliChannel, jl_convert_multi_qubit_pauli_channel)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ti", Ti, jl_convert_ti)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CV", CV, jl_convert_cv)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:StartVerbatimBox", StartVerbatimBox, jl_convert_startverbatim)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ECR", ECR, jl_convert_ecr)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CSwap", CSwap, jl_convert_cswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Ry", Ry, jl_convert_ry)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CY", CY, jl_convert_cy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CCNot", CCNot, jl_convert_ccnot)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PauliChannel", PauliChannel, jl_convert_paulichannel)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:I", I, jl_convert_i)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Unitary", Unitary, jl_convert_unitary)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Z", Z, jl_convert_z)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Si", Si, jl_convert_si)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift01", CPhaseShift01, jl_convert_cphaseshift01)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:AmplitudeDamping", AmplitudeDamping, jl_convert_amplitudedamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PSwap", PSwap, jl_convert_pswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:BitFlip", BitFlip, jl_convert_bitflip)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseDamping", PhaseDamping, jl_convert_phasedamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rz", Rz, jl_convert_rz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:GeneralizedAmplitudeDamping", GeneralizedAmplitudeDamping, jl_convert_generalizedampdamp)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseShift", PhaseShift, jl_convert_phaseshift)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:V", V, jl_convert_v)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XX", XX, jl_convert_xx)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Y", Y, jl_convert_y)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ZZ", ZZ, jl_convert_zz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Swap", Swap, jl_convert_swap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:ISwap", ISwap, jl_convert_iswap)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:H", H, jl_convert_h)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CPhaseShift", CPhaseShift, jl_convert_cphaseshift)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:PhaseFlip", PhaseFlip, jl_convert_phaseflip)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:S", S, jl_convert_s)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Depolarizing", Depolarizing, jl_convert_depo)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Rx", Rx, jl_convert_rx)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:YY", YY, jl_convert_yy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:EndVerbatimBox", EndVerbatimBox, jl_convert_endverbatim)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:T", T, jl_convert_t)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:CZ", CZ, jl_convert_cz)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:XY", XY, jl_convert_xy)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:Vi", Vi, jl_convert_vi)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Sample", Braket.IR.Sample, jl_convert_sample)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Expectation", Braket.IR.Expectation, jl_convert_expectation)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Variance", Braket.IR.Variance, jl_convert_variance)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Probability", Braket.IR.Probability, jl_convert_probability)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:StateVector", Braket.IR.StateVector, jl_convert_statevector)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:DensityMatrix", Braket.IR.DensityMatrix, jl_convert_densitymatrix)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Amplitude", Braket.IR.Amplitude, jl_convert_amplitude)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Amplitude", AbstractProgramResult, jl_convert_amplitude)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Expectation", AbstractProgramResult, jl_convert_expectation)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Probability", AbstractProgramResult, jl_convert_probability)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Sample", AbstractProgramResult, jl_convert_sample)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:StateVector", AbstractProgramResult, jl_convert_statevector)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:DensityMatrix", AbstractProgramResult, jl_convert_densitymatrix)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Variance", AbstractProgramResult, jl_convert_variance)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:AdjointGradient", AbstractProgramResult, jl_convert_adjointgradient)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:AdjointGradient", AdjointGradient, jl_convert_adjointgradient)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Amplitude", Braket.Result, jl_convert_amplitude)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Expectation", Braket.Result, jl_convert_expectation)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Probability", Braket.Result, jl_convert_probability)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Sample", Braket.Result, jl_convert_sample)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:StateVector", Braket.Result, jl_convert_statevector)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:DensityMatrix", Braket.Result, jl_convert_densitymatrix)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:Variance", Braket.Result, jl_convert_variance)
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.results:AdjointGradient", Braket.Result, jl_convert_adjointgradient)
    PythonCall.pyconvert_add_rule("braket.ir.openqasm.program_v1:Program", OpenQasmProgram, jl_convert_oqprogram)
end

function BraketSimulator.parse_program(simulator::D, program::OpenQasmProgram, shots::Int) where {D<:AbstractSimulator}
    pc      = braket[].default_simulator.openqasm.program_context.ProgramContext()
    interp  = braket[].default_simulator.openqasm.interpreter.Interpreter(pc)
    inputs  = isnothing(program.inputs) ? pybuiltins.None : pydict(program.inputs)
    py_circ = interp.build_circuit(source=pystr(program.source), inputs=inputs, is_file=endswith(program.source, ".qasm"))
    if shots > 0
        py_circ.instructions += py_circ.basis_rotation_instructions
    end
    return ir(pyconvert(Circuit, py_circ), Val(:JAQCD))
end

function simulate(
    simulator::AbstractSimulator,
    task_specs::Union{PyList{Any},NTuple{N,PyIterable},Py},
    args...;
    input::Union{PyList{Any},PyDict{Any,Any},Py}=PyDict{Any,Any}(),
    kwargs...,
) where {N}
    # handle inputs
    jl_inputs = nothing
    shots = args[end]
    stats = @timed begin
        if input isa PyDict{Any,Any}
            jl_inputs = pyconvert(Dict{String,Float64}, input)
        else
            jl_inputs = [pyconvert(Dict{String,Float64}, py_inputs) for py_inputs in input]
        end
        jl_specs = map(task_specs) do spec 
            return if pyhasattr(spec, "source")
                pyconvert(OpenQasmProgram, spec)
            else
                pyconvert(Program, spec)
            end
        end
        input = nothing
    end
    @debug "Time for conversion of specs and inputs: $(stats.time)."
    if haskey(kwargs, :measured_qubits)
        jl_measured_qubits = [pyconvert(Int, qubit) for qubit in kwargs[:measured_qubits]]
        kwargs = merge(Dict(kwargs...), Dict(:measured_qubits=>jl_measured_qubits))
    end
    PythonCall.GC.disable()
    if length(jl_specs) == 1
        result     = simulate(simulator, jl_specs[1], args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        py_result  = Py(result, task_specs[0])
    else # this is a batch! use a Braket.jl LocalSimulator to take advantage of thread migration
        local_sim   = Braket.LocalSimulator(simulator) 
        task_batch  = simulate(local_sim, jl_specs, args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        raw_results = results(task_batch)
        # now have to convert back to GateModelTaskResult from GateModelQuantumTaskResult
        processed_results = map(zip(raw_results, task_specs)) do (result, task_spec)
            header = Braket.braketSchemaHeader("braket.task_result.gate_model_task_result", "1")
            return Py(Braket.GateModelTaskResult(header, result.measurements, result.measurement_probabilities, result.result_types, result.measured_qubits, result.task_metadata, result.additional_metadata), task_spec)
        end
        py_result = pylist(processed_results)
    end
    PythonCall.GC.enable()
    return py_result
end

include("precompile.jl")
end
