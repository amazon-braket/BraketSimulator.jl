# OpenQASM 3 Braket Standard Gates
builtin_gates() = Dict{String, BuiltinGateDefinition}(
    # identity gate
    "i"=>BuiltinGateDefinition("i", String[], ["a"], (type="i", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # phase gate
    "phaseshift"=>BuiltinGateDefinition("phaseshift", ["λ"], ["a"], (type="phaseshift", arguments=InstructionArgument[:λ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # pauli X gate
    "x"=>BuiltinGateDefinition("x", String[], ["a"], (type="x", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # pauli Y gate
    "y"=>BuiltinGateDefinition("y", String[], ["a"], (type="y", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # pauli Z gate
    "z"=>BuiltinGateDefinition("z", String[], ["a"], (type="z", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # Hadamard gate
    "h"=>BuiltinGateDefinition("h", String[], ["a"], (type="h", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # S gate
    "s"=>BuiltinGateDefinition("s", String[], ["a"], (type="s", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # Si gate
    "si"=>BuiltinGateDefinition("si", String[], ["a"], (type="si", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # T gate
    "t"=>BuiltinGateDefinition("t", String[], ["a"], (type="t", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # Ti gate
    "ti"=>BuiltinGateDefinition("ti", String[], ["a"], (type="ti", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # V gate
    "v"=>BuiltinGateDefinition("v", String[], ["a"], (type="v", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # Vi gate
    "vi"=>BuiltinGateDefinition("vi", String[], ["a"], (type="vi", arguments=InstructionArgument[], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # RotX gate
    "rx"=>BuiltinGateDefinition("rx", ["θ"], ["a"], (type="rx", arguments=InstructionArgument[:θ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # RotY gate
    "ry"=>BuiltinGateDefinition("ry", ["θ"], ["a"], (type="ry", arguments=InstructionArgument[:θ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # RotZ gate
    "rz"=>BuiltinGateDefinition("rz", ["θ"], ["a"], (type="rz", arguments=InstructionArgument[:θ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # CNot gate
    "cnot"=>BuiltinGateDefinition("cnot", String[], ["a", "b"], (type="cnot", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # CY gate
    "cy"=>BuiltinGateDefinition("cy", String[], ["a", "b"], (type="cy", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # CZ gate
    "cz"=>BuiltinGateDefinition("cz", String[], ["a", "b"], (type="cz", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # CV gate
    "cv"=>BuiltinGateDefinition("cv", String[], ["a", "b"], (type="cv", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # controlled-phase
    "cphaseshift"=>BuiltinGateDefinition("cphaseshift", ["λ"], ["a", "b"], (type="cphaseshift", arguments=InstructionArgument[:λ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # controlled-phase-00
    "cphaseshift00"=>BuiltinGateDefinition("cphaseshift00", ["λ"], ["a", "b"], (type="cphaseshift00", arguments=InstructionArgument[:λ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # controlled-phase-01
    "cphaseshift01"=>BuiltinGateDefinition("cphaseshift01", ["λ"], ["a", "b"], (type="cphaseshift01", arguments=InstructionArgument[:λ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # controlled-phase-10
    "cphaseshift10"=>BuiltinGateDefinition("cphaseshift10", ["λ"], ["a", "b"], (type="cphaseshift10", arguments=InstructionArgument[:λ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # Swap gate
    "swap"=>BuiltinGateDefinition("swap", String[], ["a", "b"], (type="swap", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # ISwap gate
    "iswap"=>BuiltinGateDefinition("iswap", String[], ["a", "b"], (type="iswap", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # ISwap gate
    "pswap"=>BuiltinGateDefinition("pswap", ["θ"], ["a", "b"], (type="pswap", arguments=InstructionArgument[:θ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # controlled-swap gate
    "cswap"=>BuiltinGateDefinition("cswap", String[], ["a", "b", "c"], (type="cswap", arguments=InstructionArgument[], targets=[0, 1, 2], controls=Pair{Int,Int}[], exponent=1.0)),
    # ccnot/Toffoli gate
    "ccnot"=>BuiltinGateDefinition("ccnot", String[], ["a", "b", "c"], (type="ccnot", arguments=InstructionArgument[], targets=[0, 1, 2], controls=Pair{Int,Int}[], exponent=1.0)),
    # XX gate
    "xx"=>BuiltinGateDefinition("xx", ["θ"], ["a", "b"], (type="xx", arguments=InstructionArgument[:θ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # XY gate
    "xy"=>BuiltinGateDefinition("xy", ["θ"], ["a", "b"], (type="xy", arguments=InstructionArgument[:θ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # YY gate
    "yy"=>BuiltinGateDefinition("yy", ["θ"], ["a", "b"], (type="yy", arguments=InstructionArgument[:θ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # ZZ gate
    "zz"=>BuiltinGateDefinition("zz", ["θ"], ["a", "b"], (type="zz", arguments=InstructionArgument[:θ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # ECR gate
    "ecr"=>BuiltinGateDefinition("ecr", String[], ["a", "b"], (type="ecr", arguments=InstructionArgument[], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # MS gate
    "ms"=>BuiltinGateDefinition("ms", ["ϕ", "θ", "λ"], ["a", "b"], (type="ms", arguments=InstructionArgument[:ϕ, :θ, :λ], targets=[0, 1], controls=Pair{Int,Int}[], exponent=1.0)),
    # GPi gate
    "gpi"=>BuiltinGateDefinition("gpi", ["θ"], ["a"], (type="gpi", arguments=InstructionArgument[:θ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # GPi2 gate
    "gpi2"=>BuiltinGateDefinition("gpi2", ["θ"], ["a"], (type="gpi2", arguments=InstructionArgument[:θ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # PRx gate
    "prx"=>BuiltinGateDefinition("prx", ["θ", "ϕ"], ["a"], (type="prx", arguments=InstructionArgument[:θ, :ϕ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
    # 3-angle U gate
    "U"=>BuiltinGateDefinition("U", ["θ", "ϕ", "λ"], ["a"], (type="u", arguments=InstructionArgument[:θ, :ϕ, :λ], targets=[0], controls=Pair{Int,Int}[], exponent=1.0)),
)
