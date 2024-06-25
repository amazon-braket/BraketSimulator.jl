using BraketSimulator, PythonCall, BenchmarkTools

max_qubits(::Val{:noise}) = 5
max_qubits(::Val{:pure})  = 10

gate_operations  = pyimport("braket.default_simulator.gate_operations")
noise_operations = pyimport("braket.default_simulator.noise_operations")
local_sv         = pyimport("braket.default_simulator.state_vector_simulation")
local_dm         = pyimport("braket.default_simulator.density_matrix_simulation")
qml              = pyimport("pennylane")
np               = pyimport("numpy")
pnp              = pyimport("pennylane.numpy")
suite          = BenchmarkGroup()
suite["gates"] = BenchmarkGroup()
suite["noise"] = BenchmarkGroup()
for gate in [
    "H",
    "X",
    "Y",
    "Z",
    "V",
    "Vi",
    "T",
    "Ti",
    "S",
    "Si",
    "Rx",
    "Ry",
    "Rz",
    "GPi",
    "GPi2",
    "MS",
    "PhaseShift",
    "CNot",
    "CY",
    "CZ",
    "CV",
    "XX",
    "XY",
    "YY",
    "ZZ",
    "ECR",
    "Swap",
    "ISwap",
    "PSwap",
    "CCNot",
    "CSwap",
    "CPhaseShift",
    "CPhaseShift00",
    "CPhaseShift01",
    "CPhaseShift10",
]
    suite["gates"][gate] = BenchmarkGroup()
end
for noise in [
    "BitFlip",
    "PhaseFlip",
    "Depolarizing",
    "PauliChannel",
    "AmplitudeDamping",
    "GeneralizedAmplitudeDamping",
    "PhaseDamping",
    "TwoQubitDepolarizing",
    "TwoQubitPauliChannel",
    "TwoQubitDephasing",
    "Kraus",
]
    suite["noise"][noise] = BenchmarkGroup()
end
for n_qubits = 4:2:max_qubits(Val(:pure))
    n_amps   = 2^n_qubits
    angle    = π / 3.5
    angle2   = π / 5.2
    angle3   = π / 0.6
    for q in 0:0#n_qubits-1
        for (gate_str, gate, py_gate) in zip(
            ["H", "X", "Y", "Z", "V", "Vi", "T", "Ti", "S", "Si"],
            [BraketSimulator.H(),
             BraketSimulator.X(),
             BraketSimulator.Y(),
             BraketSimulator.Z(),
             BraketSimulator.V(),
             BraketSimulator.Vi(),
             BraketSimulator.T(),
             BraketSimulator.Ti(),
             BraketSimulator.S(),
             BraketSimulator.Si(),
            ],
            [
                [gate_operations.Hadamard([q])],
                [gate_operations.PauliX([q])],
                [gate_operations.PauliY([q])],
                [gate_operations.PauliZ([q])],
                [gate_operations.V([q])],
                [gate_operations.Vi([q])],
                [gate_operations.T([q])],
                [gate_operations.Ti([q])],
                [gate_operations.S([q])],
                [gate_operations.Si([q])],
            ],
        )
            suite["gates"][gate_str][(string(n_qubits), string(q), "Julia")] =
                @benchmarkable BraketSimulator.apply_gate!($gate, sv, $q) setup =
                    (sv = zeros(ComplexF64, $n_amps))
            suite["gates"][gate_str][(string(n_qubits), string(q), "BraketSV")] =
                @benchmarkable py_sv.evolve($py_gate) setup =
                    (py_sv = local_sv.StateVectorSimulation($n_qubits, 0, 1))
        end
        for (gate_str, gate, py_gate) in zip(
            ["Rx", "Ry", "Rz", "PhaseShift", "GPi", "GPi2"],
            [BraketSimulator.Rx(angle),
             BraketSimulator.Ry(angle),
             BraketSimulator.Rz(angle),
             BraketSimulator.PhaseShift(angle),
             BraketSimulator.GPi(angle),
             BraketSimulator.GPi2(angle),
            ],
            [
                gate_operations.RotX([q], angle),
                gate_operations.RotY([q], angle),
                gate_operations.RotZ([q], angle),
                gate_operations.PhaseShift([q], angle),
                gate_operations.GPi([q], angle),
                gate_operations.GPi2([q], angle),
            ],
        )
            suite["gates"][gate_str][(string(n_qubits), string(q), "Julia")] =
                @benchmarkable BraketSimulator.apply_gate!($gate, sv, $q) setup =
                    (sv = zeros(ComplexF64, $n_amps))
            suite["gates"][gate_str][(string(n_qubits), string(q), "BraketSV")] =
                @benchmarkable py_sv.evolve([$py_gate]) setup =
                    (py_sv = local_sv.StateVectorSimulation($n_qubits, 0, 1))
        end
        for q2 in [mod(q+1, n_qubits)] #setdiff(0:n_qubits-1, q)
            for (gate_str, gate, py_gate) in zip(
                [
                    "XX",
                    "XY",
                    "YY",
                    "ZZ",
                    "CPhaseShift",
                    "CPhaseShift00",
                    "CPhaseShift10",
                    "CPhaseShift01",
                    "PSwap",
                    "MS",
                ],
                [
                    BraketSimulator.XX(angle),
                    BraketSimulator.XY(angle),
                    BraketSimulator.YY(angle),
                    BraketSimulator.ZZ(angle),
                    BraketSimulator.CPhaseShift(angle),
                    BraketSimulator.CPhaseShift00(angle),
                    BraketSimulator.CPhaseShift10(angle),
                    BraketSimulator.CPhaseShift01(angle),
                    BraketSimulator.PSwap(angle),
                    BraketSimulator.MS(angle, angle2, angle3),
                ],
                [
                    gate_operations.XX([q, q2], angle),
                    gate_operations.XY([q, q2], angle),
                    gate_operations.YY([q, q2], angle),
                    gate_operations.ZZ([q, q2], angle),
                    gate_operations.CPhaseShift([q, q2], angle),
                    gate_operations.CPhaseShift00([q, q2], angle),
                    gate_operations.CPhaseShift10([q, q2], angle),
                    gate_operations.CPhaseShift01([q, q2], angle),
                    gate_operations.PSwap([q, q2], angle),
                    gate_operations.MS([q, q2], angle, angle2, angle3),
                ],
            )
                suite["gates"][gate_str][(string(n_qubits), string(q), string(q2), "Julia")] =
                    @benchmarkable BraketSimulator.apply_gate!($gate, sv, $q, $q2) setup =
                        (sv = zeros(ComplexF64, $n_amps))
                suite["gates"][gate_str][(string(n_qubits), string(q), string(q2), "BraketSV")] =
                    @benchmarkable py_sv.evolve([$py_gate]) setup =
                        (py_sv = local_sv.StateVectorSimulation($n_qubits, 0, 1))
            end
            for (gate_str, gate, py_gate) in zip(
                ["CNot", "CY", "CZ", "CV", "Swap", "ISwap", "ECR"],
                [BraketSimulator.CNot(),
                 BraketSimulator.CY(),
                 BraketSimulator.CZ(),
                 BraketSimulator.CV(),
                 BraketSimulator.Swap(),
                 BraketSimulator.ISwap(),
                 BraketSimulator.ECR(),
                ],
                [
                    gate_operations.CX([q, q2]),
                    gate_operations.CY([q, q2]),
                    gate_operations.CZ([q, q2]),
                    gate_operations.CV([q, q2]),
                    gate_operations.Swap([q, q2]),
                    gate_operations.ISwap([q, q2]),
                    gate_operations.ECR([q, q2]),
                ],
            )
                suite["gates"][gate_str][(string(n_qubits), string(q), string(q2), "Julia")] =
                    @benchmarkable BraketSimulator.apply_gate!($gate, sv, $q, $q2) setup =
                        (sv = zeros(ComplexF64, $n_amps))
                suite["gates"][gate_str][(string(n_qubits), string(q), string(q2), "BraketSV")] =
                    @benchmarkable py_sv.evolve([$py_gate]) setup =
                        (py_sv = local_sv.StateVectorSimulation($n_qubits, 0, 1))
            end
            for q3 in [mod(q+2, n_qubits)] #setdiff(0:n_qubits - 1, q, q2)
                for (gate_str, gate, py_gate) in zip(
                    ["CCNot", "CSwap"],
                    [BraketSimulator.CCNot(), BraketSimulator.CSwap()],
                    [gate_operations.CCNot([q, q2, 2]), gate_operations.CSwap([q, q2, 2])],
                )
                    suite["gates"][gate_str][(
                        string(n_qubits),
                        string(q),
                        string(q2),
                        string(q3),
                        "Julia",
                    )] = @benchmarkable BraketSimulator.apply_gate!($gate, sv, $q, $q2, $q3) setup =
                        (sv = zeros(ComplexF64, $n_amps))
                    suite["gates"][gate_str][(
                        string(n_qubits),
                        string(q),
                        string(q2),
                        string(q3),
                        "BraketSV",
                    )] = @benchmarkable py_sv.evolve([$py_gate]) setup =
                        (py_sv = local_sv.StateVectorSimulation($n_qubits, 0, 1))
                end
            end
        end
    end
end

for n_qubits = 2:2:max_qubits(Val(:noise))
    n_amps = 2^n_qubits
    prob = 0.1
    gamma = 0.2
    for q in 0:0#n_qubits-1
        for (noise_str, noise, py_noise) in zip(
            ["BitFlip", "PhaseFlip", "Depolarizing", "AmplitudeDamping", "PhaseDamping"],
            [
                BraketSimulator.BitFlip(prob),
                BraketSimulator.PhaseFlip(prob),
                BraketSimulator.Depolarizing(prob),
                BraketSimulator.AmplitudeDamping(prob),
                BraketSimulator.PhaseDamping(prob),
            ],
            [
                noise_operations.BitFlip([q], prob),
                noise_operations.PhaseFlip([q], prob),
                noise_operations.Depolarizing([q], prob),
                noise_operations.AmplitudeDamping([q], prob),
                noise_operations.PhaseDamping([q], prob),
            ],
        )
            suite["noise"][noise_str][(string(n_qubits), string(q), "Julia")] =
                @benchmarkable BraketSimulator.apply_noise!($noise, dm, $q) setup =
                    (dm = zeros(ComplexF64, $n_amps, $n_amps))
            suite["noise"][noise_str][(string(n_qubits), string(q), "BraketSV")] =
                @benchmarkable py_sv.evolve([$py_noise]) setup =
                    (py_sv = local_dm.DensityMatrixSimulation($n_qubits, 0))
        end
        for q2 in [mod(q+1, n_qubits)] #setdiff(0:n_qubits-1, q)
            for (noise_str, noise, py_noise) in zip(
                ["TwoQubitDepolarizing", "TwoQubitDephasing"],
                [BraketSimulator.TwoQubitDepolarizing(prob),
                 BraketSimulator.TwoQubitDephasing(prob),
                ],
                [
                    noise_operations.TwoQubitDepolarizing([q, q2], prob),
                    noise_operations.TwoQubitDephasing([q, q2], prob),
                ],
            )
                suite["noise"][noise_str][(string(n_qubits), string(q), string(q2), "Julia")] =
                    @benchmarkable BraketSimulator.apply_noise!($noise, dm, $q, $q2) setup =
                        (dm = zeros(ComplexF64, $n_amps, $n_amps))
                suite["noise"][noise_str][(string(n_qubits), string(q), string(q2), "BraketSV")] =
                    @benchmarkable py_sv.evolve([$py_noise]) setup =
                        (py_sv = local_dm.DensityMatrixSimulation($n_qubits, 0))
            end
        end
        suite["noise"]["GeneralizedAmplitudeDamping"][(string(n_qubits), string(q), "Julia")] =
            @benchmarkable BraketSimulator.apply_noise!(
                BraketSimulator.GeneralizedAmplitudeDamping($prob, $gamma),
                dm,
                $q,
            ) setup = (dm = zeros(ComplexF64, $n_amps, $n_amps))
        suite["noise"]["GeneralizedAmplitudeDamping"][(string(n_qubits), string(q), "BraketSV")] =
            @benchmarkable py_sv.evolve([
                noise_operations.GeneralizedAmplitudeDamping([$q], $prob, $gamma),
            ]) setup = (py_sv = local_dm.DensityMatrixSimulation($n_qubits, 0))
        suite["noise"]["PauliChannel"][(string(n_qubits), string(q), "Julia")] =
            @benchmarkable BraketSimulator.apply_noise!(
                PauliChannel($prob, $gamma, 0.0),
                dm,
                $q,
            ) setup = (dm = zeros(ComplexF64, $n_amps, $n_amps))
        suite["noise"]["PauliChannel"][(string(n_qubits), string(q), "BraketSV")] =
            @benchmarkable py_sv.evolve([
                noise_operations.PauliChannel(
                    targets = [$q],
                    probX = $prob,
                    probY = $gamma,
                    probZ = 0.0,
                ),
            ]) setup = (py_sv = local_dm.DensityMatrixSimulation($n_qubits, 0))
    end
end
# this is expensive! only do it if we're sure we need to regen parameters
if !isfile("small_gate_kernel_params.json")
    tune!(suite; verbose=true)
    BenchmarkTools.save("small_gate_kernel_params.json", params(suite))
end
loadparams!(suite, BenchmarkTools.load("small_gate_kernel_params.json")[1], :evals, :samples);
results = run(suite; verbose = true)
BenchmarkTools.save("small_gate_kernel_results.json", results)
