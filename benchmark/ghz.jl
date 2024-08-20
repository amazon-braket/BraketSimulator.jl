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
qiskit           = pyimport("qiskit")
qiskit_aer       = pyimport("qiskit_aer.primitives")

function ghz_circuit(qubit_count::Int, shots::Int, noise::Bool)
    ghz_circ = BraketSimulator.Circuit()
    push!(ghz_circ.instructions, BraketSimulator.Instruction(BraketSimulator.H(), 0))
    noise && push!(ghz_circ.instructions, BraketSimulator.Instruction(BraketSimulator.BitFlip(0.1), 0))
    for target_qubit = 1:qubit_count-1
        push!(ghz_circ.instructions, BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(0, target_qubit)))
        if noise
            push!(ghz_circ.instructions, BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.1), 0))
            push!(ghz_circ.instructions, BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.1), 1))
        end
    end
    shots == 0 && noise && (push!(ghz_circ.result_types, BraketSimulator.DensityMatrix()))
    shots == 0 && !noise && (push!(ghz_circ.result_types, BraketSimulator.StateVector()))
    return ghz_circ
end

function braket_sv_ghz(qubit_count::Int)
    ghz_ops = []
    push!(ghz_ops, gate_operations.Hadamard([0]))
    for target_qubit in 1:qubit_count-1
        push!(ghz_ops, gate_operations.CNot([0, target_qubit]))
    end
    return ghz_ops
end

function pl_ghz(qubit_count::Int)
    ops = [qml.Hadamard(0)]
    for target_qubit in 1:qubit_count -1
        push!(ops, qml.CNOT([0, target_qubit]))
    end
    return qml.tape.QuantumTape(ops, [])
end

function qiskit_ghz(qubit_count::Int)
    qc = qiskit.QuantumCircuit(qubit_count)
    qc.h(0)
    for target_qubit in 1:qubit_count -1
        qc.cx(0, target_qubit)
    end
    qc.measure_all()
    return qc
end

for mode in (:noise, :pure), n_qubits in 4:max_qubits(Val(mode))
    g = BenchmarkGroup()
    for shots in (0, 100)
        pl_shots = shots == 0 ? nothing : shots
        if mode == :pure
            g["Julia"]           = @benchmarkable simulate(sim, circ, $shots) setup = (sim=StateVectorSimulator($n_qubits, $shots); circ = BraketSimulator.Program(ghz_circuit($n_qubits, $shots, false)))
        elseif mode == :noise
            g["Julia"]           = @benchmarkable simulate(sim, circ, $shots) setup = (sim=DensityMatrixSimulator($n_qubits, $shots); circ = BraketSimulator.Program(ghz_circuit($n_qubits, $shots, true)))
        end
        if shots > 0 && mode == :pure && n_qubits <= 30
            # disabled for now as the existing Python simulator is very slow above 20 qubits 
            #=if n_qubits <= 27
                g["BraketSV"]    = @benchmarkable sim.execute(circ) setup = (sim=qml.device("braket.local.qubit", backend="braket_sv", wires=$n_qubits, shots=$shots); circ = pl_ghz($n_qubits))
            end=#
            g["Lightning.Qubit"] = @benchmarkable sim.execute(circ) setup = (sim=qml.device("lightning.qubit", wires=$n_qubits, shots=$pl_shots); circ = pl_ghz($n_qubits))
            g["Qiskit.Aer"]      = @benchmarkable sim.run(circ, shots=$shots).result() setup = (sim=qiskit_aer.SamplerV2(); circ = [qiskit_ghz($n_qubits)])
        end
    end
    root_fn   = "ghz_" * string(n_qubits) * "_" * string(mode) * "_"
    param_fn  = root_fn * "params.json"
    result_fn = root_fn * "results.json"
    # this is expensive! only do it if we're sure we need to regen parameters
    if !isfile(param_fn)
        tune!(g; verbose=true)
        BenchmarkTools.save(param_fn, params(g))
    end
    loadparams!(g, BenchmarkTools.load(param_fn)[1], :evals, :samples);
    results = run(g; verbose = true)
    BenchmarkTools.save(result_fn, results)
end

