if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let
        for (instructions, qubit_count, state_vector, probability_amplitudes) in [
            ([Instruction(H(), [0])], 1, [0.70710678, 0.70710678], [0.5, 0.5]),
            ([Instruction(X(), [0])], 1, [0, 1], [0, 1]),
            ([Instruction(X(), [0])], 2, [0, 0, 1, 0], [0, 0, 1, 0]),
            ([Instruction(Y(), [0])], 1, [0, im], [0, 1]),
            ([Instruction(X(), [0]), Instruction(X(), [1])], 2, [0, 0, 0, 1], [0, 0, 0, 1]),
            ([Instruction(X(), [0]), Instruction(Z(), [0])], 1, [0, -1], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(CNot(), [0, 1])],
                2,
                [0, 0, 0, 1],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(CY(), [0, 1])],
                2,
                [0, 0, 0, im],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(CZ(), [0, 1])],
                2,
                [0, 0, 1, 0],
                [0, 0, 1, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(Swap(), [0, 1])],
                2,
                [0, 1, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(ISwap(), [0, 1])],
                2,
                [0, im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(Swap(), [0, 2])],
                3,
                [0, 1, 0, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0],
            ),
            ([Instruction(X(), [0]), Instruction(S(), [0])], 1, [0, im], [0, 1]),
            ([Instruction(X(), [0]), Instruction(Si(), [0])], 1, [0, -im], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(T(), [0])],
                1,
                [0, 0.70710678 + 0.70710678 * im],
                [0, 1],
            ),
            (
                [Instruction(X(), [0]), Instruction(Ti(), [0])],
                1,
                [0, 0.70710678 - 0.70710678 * im],
                [0, 1],
            ),
            ([Instruction(V(), [0])], 1, [0.5 + 0.5 * im, 0.5 - 0.5 * im], [0.5, 0.5]),
            ([Instruction(Vi(), [0])], 1, [0.5 - 0.5 * im, 0.5 + 0.5 * im], [0.5, 0.5]),
            ([Instruction(I(), [0])], 1, [1, 0], [1, 0]),
            ([Instruction(Unitary([0 1; 1 0]), [0])], 1, [0, 1], [0, 1]),
            (
                [Instruction(X(), [0]), Instruction(PhaseShift(0.15), [0])],
                1,
                [0, 0.98877108 + 0.14943813 * im],
                [0, 1],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CPhaseShift(0.15), [0, 1]),
                ],
                2,
                [0, 0, 0, 0.98877108 + 0.14943813 * im],
                [0, 0, 0, 1],
            ),
            (
                [Instruction(CPhaseShift00(0.15), [0, 1])],
                2,
                [0.98877108 + 0.14943813 * im, 0, 0, 0],
                [1, 0, 0, 0],
            ),
            (
                [Instruction(X(), [1]), Instruction(CPhaseShift01(0.15), [0, 1])],
                2,
                [0, 0.98877108 + 0.14943813 * im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(CPhaseShift10(0.15), [0, 1])],
                2,
                [0, 0, 0.98877108 + 0.14943813 * im, 0],
                [0, 0, 1, 0],
            ),
            (
                [Instruction(Rx(0.15), [0])],
                1,
                [0.99718882, -0.07492971 * im],
                [0.99438554, 0.00561446],
            ),
            (
                [Instruction(X(), [0]), Instruction(Ry(0.15), [0])],
                1,
                [-0.07492971, 0.99718882],
                [0.00561446, 0.99438554],
            ),
            (
                [Instruction(H(), [0]), Instruction(Rz(0.15), [0])],
                1,
                [0.70511898 - 0.0529833 * im, 0.70511898 + 0.0529833 * im],
                [0.5, 0.5],
            ),
            (
                [Instruction(X(), [0]), Instruction(PSwap(0.15), [0, 1])],
                2,
                [0, 0.98877108 + 0.14943813 * im, 0, 0],
                [0, 1, 0, 0],
            ),
            (
                [Instruction(X(), [0]), Instruction(XY(0.15), [0, 1])],
                2,
                [0, 0.07492971 * im, 0.99718882, 0],
                [0, 0.00561446, 0.99438554, 0],
            ),
            (
                [Instruction(XX(0.3), [0, 1])],
                2,
                [0.98877108, 0, 0, -0.14943813 * im],
                [0.97766824, 0, 0, 0.02233176],
            ),
            (
                [Instruction(YY(0.3), [0, 1])],
                2,
                [0.98877108, 0, 0, 0.14943813 * im],
                [0.97766824, 0, 0, 0.02233176],
            ),
            (
                [Instruction(ZZ(0.15), [0, 1])],
                2,
                [0.99718882 - 0.07492971 * im, 0, 0, 0],
                [1, 0, 0, 0],
            ),
            (
                [Instruction(BraketSimulator.MultiRZ(0.15), [0, 1, 2])],
                3,
                [0.99718882 - 0.07492971 * im, 0, 0, 0, 0, 0, 0, 0],
                [1, 0, 0, 0, 0, 0, 0, 0],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CCNot(), [0, 1, 2]),
                ],
                3,
                [0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 0, 0, 0, 0, 0, 1],
            ),
            (
                [
                    Instruction(X(), [0]),
                    Instruction(X(), [1]),
                    Instruction(CSwap(), [0, 1, 2]),
                ],
                3,
                [0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0],
            ),
        ]
            simulation = StateVectorSimulator(qubit_count, 0)
            simulation = evolve!(simulation, instructions)
        end
        for (obs, equivalent_gates, qubit_count) in [
            ([(Braket.Observables.X(), [0])], [Instruction(H(), [0])], 1),
            ([(Braket.Observables.Z(), [0])], Instruction[], 1),
            ([(Braket.Observables.I(), [0])], Instruction[], 1),
            (
                [
                    (Braket.Observables.X(), [0]),
                    (Braket.Observables.Z(), [3]),
                    (Braket.Observables.H(), [2]),
                ],
                [Instruction(H(), [0]), Instruction(Ry(-π / 4), [2])],
                5,
            ),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.X(),
                        Braket.Observables.Z(),
                        Braket.Observables.H(),
                        Braket.Observables.I(),
                    ]),
                    (0, 3, 2, 1),
                )],
                [Instruction(H(), [0]), Instruction(Ry(-π / 4), [2])],
                5,
            ),
            (
                [(Braket.Observables.X(), [0, 1])],
                [Instruction(H(), [0]), Instruction(H(), [1])],
                2,
            ),
            ([(Braket.Observables.Z(), [0, 1])], Instruction[], 2),
            ([(Braket.Observables.I(), [0, 1])], Instruction[], 2),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.I(),
                        Braket.Observables.Z(),
                    ]),
                    (2, 0),
                )],
                Instruction[],
                3,
            ),
            (
                [(
                    Braket.Observables.TensorProduct([
                        Braket.Observables.X(),
                        Braket.Observables.Z(),
                    ]),
                    (2, 0),
                )],
                [Instruction(H(), [2])],
                3,
            ),
        ]
            sim_observables = StateVectorSimulator(qubit_count, 0)
            sim_observables = BraketSimulator.apply_observables!(sim_observables, obs)
            sim_gates = StateVectorSimulator(qubit_count, 0)
            sim_gates = BraketSimulator.evolve!(sim_gates, equivalent_gates)
            @assert BraketSimulator.state_with_observables(sim_observables) ≈
                    sim_gates.state_vector
        end

        function qft_circuit_operations(qubit_count::Int)
            qft_ops = Instruction[]
            for target_qubit = 0:qubit_count-1
                angle = π / 2
                push!(qft_ops, Instruction(H(), [target_qubit]))
                for control_qubit = target_qubit+1:qubit_count-1
                    push!(
                        qft_ops,
                        Instruction(CPhaseShift(angle), [control_qubit, target_qubit]),
                    )
                    angle /= 2
                end
            end
            return qft_ops
        end

        qubit_count = 16
        simulation = StateVectorSimulator(qubit_count, 0)
        operations = qft_circuit_operations(qubit_count)
        simulation = BraketSimulator.evolve!(simulation, operations)
        @assert BraketSimulator.probabilities(simulation) ≈
                fill(1.0 / (2^qubit_count), 2^qubit_count)
        simulation = StateVectorSimulator(2, 10000)
        simulation = BraketSimulator.evolve!(
            simulation,
            [Instruction(H(), [0]), Instruction(CNot(), [0, 1])],
        )
        samples = BraketSimulator.samples(simulation)

        for N = 1:16
            pe = Braket.PauliEigenvalues(Val(N))
            for i = 1:2^N
                pe[i]
            end
        end
        # optional PythonCall stuff
        neededPackages = [:PythonCall] 
        using Pkg
        for neededpackage in neededPackages
            if String(neededpackage) in keys(Pkg.project().dependencies)
                @eval using $neededpackage
                n_qubits = 10
                c = braket[].circuits.circuit.Circuit()
                c.h(0)
                # ghz
                for q in 1:n_qubits-1
                    c.cnot(0, q)
                end
                c.state_vector()
                svs = StateVectorSimulator(0, 0)
                svs([c], shots=0)
            end
        end
    end
end
