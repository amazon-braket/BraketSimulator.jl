suite["ghz"] = BenchmarkGroup()

function ghz_circuit(qubit_count::Int)
    ghz_circ = Circuit()
    H(ghz_circ, 0)
    for target_qubit = 1:qubit_count-1
        CNot(ghz_circ, 0, target_qubit)
    end
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
    return qml.tape.QuantumTape(ops, [qml.sample()])
end

for mode in (:noise, :pure), n_qubits in 4:max_qubits(Val(mode)), shots in (0, 100)
    key = (string(n_qubits), string(shots), string(mode))
    suite["ghz"][key] = BenchmarkGroup()
    g = suite["ghz"][key]
    if mode == :pure
        g["Julia"]           = @benchmarkable sim(circ, shots=shots) setup = (sim=LocalSimulator("braket_sv"); circ = ghz_circuit($n_qubits))
        g["Lightning.Qubit"] = @benchmarkable sim.execute(circ) setup = (sim=qml.device("lightning.qubit", wires=$n_qubits, shots=shots); circ = pl_ghz($n_qubits))
    elseif mode == :noise
        g["Julia"]           = @benchmarkable sim(circ, shots=shots) setup = (sim=LocalSimulator("braket_sv"); circ = Braket.apply_noise(nm, ghz_circuit($n_qubits)))
    end
    if shots > 0 && n_qubits <= 30
        g["Qiskit.Aer"]      = @benchmarkable sim.backend.run(circ, shots=shots).result()  setup = (sim=qml.device("qiskit.aer", backend="aer_simulator_statevector", wires=$n_qubits, shots=shots, statevector_parallel_threshold=8); circ = sim.compile_circuits(pylist([pl_ghz($n_qubits)])))
    end
end
