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
    PythonCall.pyconvert_add_rule("braket.ir.jaqcd.instructions:EndVerbatimBox", EndVerbatimBox, jl_convert_endverbatim)
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

function simulate(
    d::AbstractSimulator,
    task_specs::Union{PyList{Any},NTuple{N,PyIterable},Py},
    args...;
    input::Union{PyList{Any},PyDict{Any,Any},Py}=PyDict{Any,Any}(),
    kwargs...,
) where {N}
    # handle inputs
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
            jl_prog = if pyhasattr(spec, "source")
                pyconvert(OpenQasmProgram, spec)
            else
                pyconvert(Program, spec)
            end
            push!(jl_specs, jl_prog)
            s_ix += 1
        end
        input = nothing
    end
    @debug "Time for conversion of specs and inputs: $(stats.time)."
    if haskey(kwargs, :measured_qubits)
        jl_mqs = [pyconvert(Int, q) for q in kwargs[:measured_qubits]]
        kwargs = merge(Dict(kwargs...), Dict(:measured_qubits=>jl_mqs))
    end

    PythonCall.GC.disable()
    if length(jl_specs) == 1
        r = simulate(d, jl_specs[1], args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        py_r = Py(r, task_specs[1])
    else
        r = simulate(d, jl_specs, args[1:end-1]...; inputs = jl_inputs, shots=shots, kwargs...)
        py_r = Py(r, task_specs)
    end
    PythonCall.GC.enable()
    return py_r
end

py_obs(o::String) = pylist([pystr(o)])
function py_obs(obs::Vector)
    raw_obs = map(obs) do o
        o isa String ? pystr(o) : pylist(pylist(pylist(o__) for o__ in o_) for o_ in o)
    end
    return pylist(raw_obs)
end
function Py(r::Braket.IR.Sample)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.Sample(targets=py_targets, observable=py_obs(r.observable), type=pystr("sample"))
end

function Py(r::Braket.IR.Expectation)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.Expectation(targets=py_targets, observable=py_obs(r.observable), type=pystr("expectation"))
end

function Py(r::Braket.IR.Variance)
    py_targets = isnothing(r.targets) ? PythonCall.pybuiltins.None : pylist(r.targets)
    return braket[].ir.jaqcd.results.Variance(observable=py_obs(r.observable), targets=py_targets, type=pystr("variance"))
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

function Py(r::GateModelTaskResult, act)
    py_measurements  = isempty(r.measurements) ? PythonCall.pybuiltins.None : pylist(pylist(meas) for meas in r.measurements)
    py_probabilities = !isnothing(r.measurementProbabilities) ? pydict(Dict(pystr(k)=>v for (k,v) in r.measurementProbabilities)) : PythonCall.pybuiltins.None
    py_qubits        = !isnothing(r.measuredQubits) ? pylist(r.measuredQubits) : PythonCall.pybuiltins.None
    py_results       = pylist(Py(rtv) for rtv in r.resultTypes)
    py_task_mtd      = braket[].task_result.task_metadata_v1.TaskMetadata(id=pystr(r.taskMetadata.id), shots=Py(r.taskMetadata.shots), deviceId=pystr(r.taskMetadata.deviceId))
    py_addl_mtd      = braket[].task_result.additional_metadata.AdditionalMetadata(action=braket[].ir.openqasm.program_v1.Program(source=""))
    return braket[].task_result.GateModelTaskResult(measurements=py_measurements, measurementProbabilities=py_probabilities, resultTypes=py_results, measuredQubits=py_qubits, taskMetadata=py_task_mtd, additionalMetadata=py_addl_mtd)
end

include("precompile.jl")
end
