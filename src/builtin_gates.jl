# OpenQASM 3 Braket Standard Gates

# not in this file: unitary
# #pragma braket unitary(<matrix>) target

const builtin_gates = Dict{String, GateDefinition}(
    # identity gate
    "i"=>GateDefinition("i", String[], ["a"], [Instruction(Braket.I(), 0)]),
    # phase gate
    "phaseshift"=>GateDefinition("phaseshift", ["λ"], ["a"], [Instruction(PhaseShift(FreeParameter(:λ)), 0)]),
    # pauli X gate
    "x"=>GateDefinition("x", String[], ["a"], [Instruction(X(), 0)]),
    # pauli Y gate
    "y"=>GateDefinition("y", String[], ["a"], [Instruction(Y(), 0)]),
    # pauli Z gate
    "z"=>GateDefinition("z", String[], ["a"], [Instruction(Z(), 0)]),
    # Hadamard gate
    "h"=>GateDefinition("h", String[], ["a"], [Instruction(H(), 0)]),
    # S gate
    "s"=>GateDefinition("s", String[], ["a"], [Instruction(S(), 0)]),
    # Si gate
    "si"=>GateDefinition("si", String[], ["a"], [Instruction(Si(), 0)]),
    # T gate
    "t"=>GateDefinition("t", String[], ["a"], [Instruction(T(), 0)]),
    # Ti gate
    "ti"=>GateDefinition("ti", String[], ["a"], [Instruction(Ti(), 0)]),
    # V gate
    "v"=>GateDefinition("v", String[], ["a"], [Instruction(V(), 0)]),
    # Vi gate
    "vi"=>GateDefinition("vi", String[], ["a"], [Instruction(Vi(), 0)]),
    # RotX gate
    "rx"=>GateDefinition("rx", ["θ"], ["a"], [Instruction(Rx(FreeParameter(:θ)), 0)]),
    # RotY gate
    "ry"=>GateDefinition("ry", ["θ"], ["a"], [Instruction(Ry(FreeParameter(:θ)), 0)]),
    # RotZ gate
    "rz"=>GateDefinition("rz", ["θ"], ["a"], [Instruction(Rz(FreeParameter(:θ)), 0)]),
    # CNot gate
    "cnot"=>GateDefinition("cnot", String[], ["a", "b"], [Instruction(CNot(), [0, 1])]),
    # CY gate
    "cy"=>GateDefinition("cy", String[], ["a", "b"], [Instruction(CY(), [0, 1])]),
    # CZ gate
    "cz"=>GateDefinition("cz", String[], ["a", "b"], [Instruction(CZ(), [0, 1])]),
    # CV gate
    "cv"=>GateDefinition("cv", String[], ["a", "b"], [Instruction(CV(), [0, 1])]),
    # controlled-phase
    "cphaseshift"=>GateDefinition("cphaseshift", ["λ"], ["a", "b"], [Instruction(CPhaseShift(FreeParameter(:λ)), [0, 1])]),
    # controlled-phase-00
    "cphaseshift00"=>GateDefinition("cphaseshift00", ["λ"], ["a", "b"], [Instruction(CPhaseShift00(FreeParameter(:λ)), [0, 1])]),
    # controlled-phase-01
    "cphaseshift01"=>GateDefinition("cphaseshift01", ["λ"], ["a", "b"], [Instruction(CPhaseShift01(FreeParameter(:λ)), [0, 1])]),
    # controlled-phase-10
    "cphaseshift10"=>GateDefinition("cphaseshift10", ["λ"], ["a", "b"], [Instruction(CPhaseShift10(FreeParameter(:λ)), [0, 1])]),
    # Swap gate
    "swap"=>GateDefinition("swap", String[], ["a", "b"], [Instruction(Swap(), [0, 1])]),
    # ISwap gate
    "iswap"=>GateDefinition("iswap", String[], ["a", "b"], [Instruction(ISwap(), [0, 1])]),
    # ISwap gate
    "pswap"=>GateDefinition("pswap", ["θ"], ["a", "b"], [Instruction(PSwap(FreeParameter(:θ)), [0, 1])]),
    # controlled-swap gate
    "cswap"=>GateDefinition("cswap", String[], ["a", "b", "c"], [Instruction(CSwap(), [0, 1, 2])]),
    # ccnot/Toffoli gate
    "ccnot"=>GateDefinition("ccnot", String[], ["a", "b", "c"], [Instruction(CCNot(), [0, 1, 2])]),
    # XX gate
    "xx"=>GateDefinition("xx", ["θ"], ["a", "b"], [Instruction(XX(FreeParameter(:θ)), [0, 1])]),
    # XY gate
    "xy"=>GateDefinition("xy", ["θ"], ["a", "b"], [Instruction(XY(FreeParameter(:θ)), [0, 1])]),
    # YY gate
    "yy"=>GateDefinition("yy", ["θ"], ["a", "b"], [Instruction(YY(FreeParameter(:θ)), [0, 1])]),
    # ZZ gate
    "zz"=>GateDefinition("zz", ["θ"], ["a", "b"], [Instruction(ZZ(FreeParameter(:θ)), [0, 1])]),
    # ECR gate
    "ecr"=>GateDefinition("ecr", String[], ["a", "b"], [Instruction(ECR(), [0, 1])]),
    # MS gate
    "ms"=>GateDefinition("ms", ["ϕ", "θ", "λ"], ["a", "b"], [Instruction(MS(FreeParameter(:ϕ), FreeParameter(:θ), FreeParameter(:λ)), [0, 1])]),
    # GPi gate
    "gpi"=>GateDefinition("gpi", ["θ"], ["a"], [Instruction(GPi(FreeParameter(:θ)), 0)]),
    # GPi2 gate
    "gpi2"=>GateDefinition("gpi2", ["θ"], ["a"], [Instruction(GPi2(FreeParameter(:θ)), 0)]),
    # PRx gate
    "prx"=>GateDefinition("prx", ["θ", "ϕ"], ["a"], [Instruction(PRx(FreeParameter(:θ), FreeParameter(:ϕ)), 0)]),
    # 3-angle U gate
    "U"=>GateDefinition("U", ["θ", "ϕ", "λ"], ["a"], [Instruction(U(FreeParameter(:θ), FreeParameter(:ϕ), FreeParameter(:λ)), 0)]),
)

const noise_types = Dict{String, Type}(
                                       "bit_flip"=>BitFlip,
                                       "phase_flip"=>PhaseFlip,
                                       "pauli_channel"=>PauliChannel,
                                       "depolarizing"=>Depolarizing,
                                       "two_qubit_depolarizing"=>TwoQubitDepolarizing,
                                       "two_qubit_dephasing"=>TwoQubitDephasing,
                                       "amplitude_damping"=>AmplitudeDamping,
                                       "generalized_amplitude_damping"=>GeneralizedAmplitudeDamping,
                                       "phase_damping"=>PhaseDamping,
                                       "kraus"=>Kraus,
                                      )
