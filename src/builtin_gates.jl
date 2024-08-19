# OpenQASM 3 Braket Standard Gates
builtin_gates() = Dict{String, GateDefinition}(
    # identity gate
    "i"=>GateDefinition("i", String[], ["a"], Instruction(BraketSimulator.I(), 0)),
    # phase gate
    "phaseshift"=>GateDefinition("phaseshift", ["λ"], ["a"], Instruction(BraketSimulator.PhaseShift(BraketSimulator.FreeParameter(:λ)), 0)),
    # pauli X gate
    "x"=>GateDefinition("x", String[], ["a"], Instruction(BraketSimulator.X(), 0)),
    # pauli Y gate
    "y"=>GateDefinition("y", String[], ["a"], Instruction(BraketSimulator.Y(), 0)),
    # pauli Z gate
    "z"=>GateDefinition("z", String[], ["a"], Instruction(BraketSimulator.Z(), 0)),
    # Hadamard gate
    "h"=>GateDefinition("h", String[], ["a"], Instruction(BraketSimulator.H(), 0)),
    # S gate
    "s"=>GateDefinition("s", String[], ["a"], Instruction(BraketSimulator.S(), 0)),
    # Si gate
    "si"=>GateDefinition("si", String[], ["a"], Instruction(BraketSimulator.Si(), 0)),
    # T gate
    "t"=>GateDefinition("t", String[], ["a"], Instruction(BraketSimulator.T(), 0)),
    # Ti gate
    "ti"=>GateDefinition("ti", String[], ["a"], Instruction(BraketSimulator.Ti(), 0)),
    # V gate
    "v"=>GateDefinition("v", String[], ["a"], Instruction(BraketSimulator.V(), 0)),
    # Vi gate
    "vi"=>GateDefinition("vi", String[], ["a"], Instruction(BraketSimulator.Vi(), 0)),
    # RotX gate
    "rx"=>GateDefinition("rx", ["θ"], ["a"], Instruction(BraketSimulator.Rx(BraketSimulator.FreeParameter(:θ)), 0)),
    # RotY gate
    "ry"=>GateDefinition("ry", ["θ"], ["a"], Instruction(BraketSimulator.Ry(BraketSimulator.FreeParameter(:θ)), 0)),
    # RotZ gate
    "rz"=>GateDefinition("rz", ["θ"], ["a"], Instruction(BraketSimulator.Rz(BraketSimulator.FreeParameter(:θ)), 0)),
    # CNot gate
    "cnot"=>GateDefinition("cnot", String[], ["a", "b"], Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(0, 1))),
    # CY gate
    "cy"=>GateDefinition("cy", String[], ["a", "b"], Instruction(BraketSimulator.CY(), BraketSimulator.QubitSet(0, 1))),
    # CZ gate
    "cz"=>GateDefinition("cz", String[], ["a", "b"], Instruction(BraketSimulator.CZ(), BraketSimulator.QubitSet(0, 1))),
    # CV gate
    "cv"=>GateDefinition("cv", String[], ["a", "b"], Instruction(BraketSimulator.CV(), BraketSimulator.QubitSet(0, 1))),
    # controlled-phase
    "cphaseshift"=>GateDefinition("cphaseshift", ["λ"], ["a", "b"], Instruction(BraketSimulator.CPhaseShift(BraketSimulator.FreeParameter(:λ)), BraketSimulator.QubitSet(0, 1))),
    # controlled-phase-00
    "cphaseshift00"=>GateDefinition("cphaseshift00", ["λ"], ["a", "b"], Instruction(BraketSimulator.CPhaseShift00(BraketSimulator.FreeParameter(:λ)), BraketSimulator.QubitSet(0, 1))),
    # controlled-phase-01
    "cphaseshift01"=>GateDefinition("cphaseshift01", ["λ"], ["a", "b"], Instruction(BraketSimulator.CPhaseShift01(BraketSimulator.FreeParameter(:λ)), BraketSimulator.QubitSet(0, 1))),
    # controlled-phase-10
    "cphaseshift10"=>GateDefinition("cphaseshift10", ["λ"], ["a", "b"], Instruction(BraketSimulator.CPhaseShift10(BraketSimulator.FreeParameter(:λ)), BraketSimulator.QubitSet(0, 1))),
    # Swap gate
    "swap"=>GateDefinition("swap", String[], ["a", "b"], Instruction(BraketSimulator.Swap(), BraketSimulator.QubitSet(0, 1))),
    # ISwap gate
    "iswap"=>GateDefinition("iswap", String[], ["a", "b"], Instruction(BraketSimulator.ISwap(), BraketSimulator.QubitSet(0, 1))),
    # ISwap gate
    "pswap"=>GateDefinition("pswap", ["θ"], ["a", "b"], Instruction(BraketSimulator.PSwap(BraketSimulator.FreeParameter(:θ)), BraketSimulator.QubitSet(0, 1))),
    # controlled-swap gate
    "cswap"=>GateDefinition("cswap", String[], ["a", "b", "c"], Instruction(BraketSimulator.CSwap(), BraketSimulator.QubitSet(0, 1, 2))),
    # ccnot/Toffoli gate
    "ccnot"=>GateDefinition("ccnot", String[], ["a", "b", "c"], Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(0, 1, 2))),
    # XX gate
    "xx"=>GateDefinition("xx", ["θ"], ["a", "b"], Instruction(BraketSimulator.XX(BraketSimulator.FreeParameter(:θ)), BraketSimulator.QubitSet(0, 1))),
    # XY gate
    "xy"=>GateDefinition("xy", ["θ"], ["a", "b"], Instruction(BraketSimulator.XY(BraketSimulator.FreeParameter(:θ)), BraketSimulator.QubitSet(0, 1))),
    # YY gate
    "yy"=>GateDefinition("yy", ["θ"], ["a", "b"], Instruction(BraketSimulator.YY(BraketSimulator.FreeParameter(:θ)), BraketSimulator.QubitSet(0, 1))),
    # ZZ gate
    "zz"=>GateDefinition("zz", ["θ"], ["a", "b"], Instruction(BraketSimulator.ZZ(BraketSimulator.FreeParameter(:θ)), BraketSimulator.QubitSet(0, 1))),
    # ECR gate
    "ecr"=>GateDefinition("ecr", String[], ["a", "b"], Instruction(BraketSimulator.ECR(), BraketSimulator.QubitSet(0, 1))),
    # MS gate
    "ms"=>GateDefinition("ms", ["ϕ", "θ", "λ"], ["a", "b"], Instruction(BraketSimulator.MS(BraketSimulator.FreeParameter(:ϕ), BraketSimulator.FreeParameter(:θ), BraketSimulator.FreeParameter(:λ)), BraketSimulator.QubitSet(0, 1))),
    # GPi gate
    "gpi"=>GateDefinition("gpi", ["θ"], ["a"], Instruction(BraketSimulator.GPi(BraketSimulator.FreeParameter(:θ)), 0)),
    # GPi2 gate
    "gpi2"=>GateDefinition("gpi2", ["θ"], ["a"], Instruction(BraketSimulator.GPi2(BraketSimulator.FreeParameter(:θ)), 0)),
    # PRx gate
    "prx"=>GateDefinition("prx", ["θ", "ϕ"], ["a"], Instruction(BraketSimulator.PRx(BraketSimulator.FreeParameter(:θ), BraketSimulator.FreeParameter(:ϕ)), 0)),
    # 3-angle U gate
    "U"=>GateDefinition("U", ["θ", "ϕ", "λ"], ["a"], Instruction(BraketSimulator.U(BraketSimulator.FreeParameter(:θ), BraketSimulator.FreeParameter(:ϕ), BraketSimulator.FreeParameter(:λ)), 0)),
)

const noise_types = Dict{String, Type}(
                                       "bit_flip"=>BraketSimulator.BitFlip,
                                       "phase_flip"=>BraketSimulator.PhaseFlip,
                                       "pauli_channel"=>BraketSimulator.PauliChannel,
                                       "depolarizing"=>BraketSimulator.Depolarizing,
                                       "two_qubit_depolarizing"=>BraketSimulator.TwoQubitDepolarizing,
                                       "two_qubit_dephasing"=>BraketSimulator.TwoQubitDephasing,
                                       "amplitude_damping"=>BraketSimulator.AmplitudeDamping,
                                       "generalized_amplitude_damping"=>BraketSimulator.GeneralizedAmplitudeDamping,
                                       "phase_damping"=>BraketSimulator.PhaseDamping,
                                       "kraus"=>BraketSimulator.Kraus,
                                      )
