using BraketSimulator, PythonCall, BenchmarkTools

max_qubits(::Val{:noise}) = 10
max_qubits(::Val{:pure})  = 20

gate_operations  = pyimport("braket.default_simulator.gate_operations")
noise_operations = pyimport("braket.default_simulator.noise_operations")
local_sv         = pyimport("braket.default_simulator.state_vector_simulation")
local_dm         = pyimport("braket.default_simulator.density_matrix_simulation")
qml              = pyimport("pennylane")
np               = pyimport("numpy")
pnp              = pyimport("pennylane.numpy")
qiskit           = pyimport("qiskit")
qiskit_aer       = pyimport("qiskit_aer.primitives")

function qft_circuit(qubit_count::Int, noise::Bool)
    qft_circ = BraketSimulator.Circuit()
    for target_qubit = 0:qubit_count-1
        angle = π / 2
        push!(qft_circ.instructions, BraketSimulator.Instruction(BraketSimulator.H(), target_qubit))
        noise && push!(qft_circ.instructions, BraketSimulator.Instruction(BraketSimulator.BitFlip(0.1), target_qubit))
        for control_qubit = target_qubit+1:qubit_count-1
            push!(qft_circ.instructions, BraketSimulator.Instruction(BraketSimulator.CPhaseShift(angle), [control_qubit, target_qubit]))
            noise && push!(qft_circ.instructions, BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.1), target_qubit))
            noise && push!(qft_circ.instructions, BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.1), control_qubit))
            angle /= 2
        end
    end
    push!(qft_circ.result_types, BraketSimulator.Probability())
    return qft_circ
end

function braket_sv_qft(qubit_count::Int, noise::Bool)
    qft_ops = []
    for target_qubit in 0:qubit_count-1
        angle = π / 2
        push!(qft_ops, gate_operations.Hadamard([target_qubit]))
        noise && push!(qft_ops, noise_operations.BitFlip(0.1, [target_qubit]))
        for control_qubit in target_qubit + 1:qubit_count-1
            push!(qft_ops, gate_operations.CPhaseShift([control_qubit, target_qubit], angle))
            noise && push!(qft_ops, noise_operations.PhaseFlip(0.1, [control_qubit]))
            noise && push!(qft_ops, noise_operations.PhaseFlip(0.1, [target_qubit]))
            angle /= 2
        end
    end
    return qft_ops
end
    
function pl_qft(qubit_count::Int)
    qft_ops = [qml.PauliX(0)]
    push!(qft_ops, qml.QFT(wires=collect(0:qubit_count-1)))
    tape = qml.tape.QuantumTape(qft_ops, [qml.expval(qml.operation.Tensor(qml.PauliZ(0), qml.PauliX(1)))])
    return tape.expand(depth=5)
end

function qiskit_qft(qubit_count::Int)
    qc = qiskit.QuantumCircuit(qubit_count)
    for target_qubit in 0:qubit_count -1
        angle = π / 2
        qc.h(target_qubit)
        for control_qubit in target_qubit + 1:qubit_count-1
            qc.cp(angle, control_qubit, target_qubit)
            angle /= 2
        end
    end
    qc.measure_all()
    return qc
end


for mode in (:noise, :pure), n_qubits in 4:max_qubits(Val(mode))
    g = BenchmarkGroup()
    for shots in (0, 100)
        if mode == :pure
            g["Julia-$shots"]           = @benchmarkable simulate(sim, circ, $shots) setup = (sim=StateVectorSimulator($n_qubits, $shots); circ = BraketSimulator.Program(qft_circuit($n_qubits, false)))
        elseif mode == :noise
            g["Julia-$shots"]           = @benchmarkable simulate(sim, circ, $shots) setup = (sim=DensityMatrixSimulator($n_qubits, $shots); circ = BraketSimulator.Program(qft_circuit($n_qubits, true)))
        end
        if shots > 0 && n_qubits <= 30 && mode == :pure
            if n_qubits <= 27
                g["BraketSV-$shots"]        = @benchmarkable sim.execute(circ) setup = (sim=qml.device("braket.local.qubit", backend="braket_sv", wires=$n_qubits, shots=$shots); circ = pl_qft($n_qubits))
            end
            g["Lightning.Qubit-$shots"] = @benchmarkable sim.execute(circ) setup = (sim=qml.device("lightning.qubit", wires=$n_qubits, shots=$shots); circ = pl_qft($n_qubits))
            g["Qiskit.Aer-$shots"]      = @benchmarkable sim.run(circ, shots=$shots).result() setup = (sim=qiskit_aer.SamplerV2(); circ = [qiskit_qft($n_qubits)])
        end
    end
    root_fn   = "qft_" * string(n_qubits) * "_" * string(mode) * "_"
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
