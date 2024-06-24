using Test, Braket, BraketSimulator, DataStructures

import Braket: I, Instruction

LARGE_TESTS = get(ENV, "BRAKET_SV_LARGE_TESTS", false)

@testset "State vector simulator" begin
    @testset "Simple circuits nq: $qubit_count $instructions" for (
        instructions,
        qubit_count,
        state_vector,
        probability_amplitudes,
    ) in [
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
            [Instruction(H(), [0]), Instruction(GPi(0.15), [0])],
            1,
            [0.69916673 - 0.10566872im, 0.69916673 + 0.10566872im],
            [0.5, 0.5],
        ),
        (
            [Instruction(H(), [0]), Instruction(GPi2(0.15), [0])],
            1,
            [0.42528093 - 0.49438554im, 0.57471907 - 0.49438554im],
            [0.42528093, 0.57471907],
        ),
        (
            [Instruction(MS(π/2, -π/4, 0.3), [0, 1])],
            2,
            [0.98877108, 0, 0, 0.10566872 - 0.10566872im],
            [0.97766824, 0, 0, 0.02233176],
        ),
        (
            [Instruction(H(), [0]), Instruction(H(), [1]), Instruction(ECR(), [0, 1])],
            2,
            [
                0.35355339 + 0.35355339im,
                0.35355339 + 0.35355339im,
                0.35355339 - 0.35355339im,
                0.35355339 - 0.35355339im,
            ],
            [0.25, 0.25, 0.25, 0.25],
        ),
        (
            [Instruction(H(), [0]), Instruction(U(0.15, 0.4, 0.7), [0])],
            1,
            [0.66459511 - 0.03413278im, 0.36864009 + 0.64903989im],
            [0.44285171, 0.55714829],
        ),
        (
            [Instruction(H(), [0]), Instruction(MultiQubitPhaseShift{1}(0.15), [0])],
            1,
            [0.69916673 + 0.10566872im, 0.69916673 + 0.10566872im],
            [0.5, 0.5],
        ),
        (
            [Instruction(H(), [0]), Instruction(PRx(0.15, 0.4), [0])],
            1,
            [0.6844863 - 0.04880085im, 0.72575165 - 0.04880085im],
            [0.47090303, 0.52909697],
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
        @testset "Simulator $sim" for sim in (StateVectorSimulator, DensityMatrixSimulator) 
            simulation = sim(qubit_count, 0)
            simulation = evolve!(simulation, instructions)
            if sim == StateVectorSimulator
                @test state_vector ≈ collect(BraketSimulator.state_vector(simulation))
            end
            @test probability_amplitudes ≈
                  collect(BraketSimulator.probabilities(simulation))
        end
    end
    @testset "Apply observables $obs" for (obs, equivalent_gates, qubit_count) in [
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
        @test BraketSimulator.state_with_observables(sim_observables) ≈
              sim_gates.state_vector
    end
    @testset "Apply observables fails at second call" begin
        simulation = StateVectorSimulator(4, 0)
        simulation = BraketSimulator.apply_observables!(
            simulation,
            [(Braket.Observables.X(), [0])],
        )
        @test_throws ErrorException BraketSimulator.apply_observables!(
            simulation,
            [(Braket.Observables.X(), [0])],
        )
    end
    @testset "state_with_observables fails before any observables are applied" begin
        simulation = StateVectorSimulator(4, 0)
        @test_throws ErrorException BraketSimulator.state_with_observables(simulation)
    end
    @testset "QFT simulation" begin
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
        max_qc = LARGE_TESTS ? 32 : 20
        @testset "Qubit count $qubit_count" for qubit_count in 4:max_qc
            simulation = StateVectorSimulator(qubit_count, 0)
            operations = qft_circuit_operations(qubit_count)
            simulation = BraketSimulator.evolve!(simulation, operations)
            @test collect(BraketSimulator.probabilities(simulation)) ≈
                  fill(1.0 / (2^qubit_count), 2^qubit_count)
        end
    end
    @testset "samples" begin
        simulation = StateVectorSimulator(2, 10000)
        simulation = BraketSimulator.evolve!(
            simulation,
            [Instruction(H(), [0]), Instruction(CNot(), [0, 1])],
        )
        samples = counter(BraketSimulator.samples(simulation))

        @test qubit_count(simulation) == 2
        @test collect(keys(samples)) == [0, 3]
        @test 0.4 < samples[0] / (samples[0] + samples[3]) < 0.6
        @test 0.4 < samples[3] / (samples[0] + samples[3]) < 0.6
        @test samples[0] + samples[3] == 10000
    end
    @testset "batch" begin
        function make_ghz(num_qubits)
            ghz = Circuit()
            ghz(H, 0)
            for ii in 0:num_qubits-2
                ghz(CNot, ii, ii+1)
            end
            return ir(ghz)
        end
        num_qubits = 5
        n_circuits = 100
        shots   = 1000
        jl_ghz  = [make_ghz(num_qubits) for ix in 1:n_circuits]
        jl_sim  = StateVectorSimulator(num_qubits, 0);
        results = simulate(jl_sim, jl_ghz, shots)
        for r in results
            @test 400 < count(m->m == fill(0, num_qubits), r.measurements) < 600
            @test 400 < count(m->m == fill(1, num_qubits), r.measurements) < 600
        end
    end
    @testset "similar, copy and copyto!" begin
        qubit_count = 10
        orig = StateVectorSimulator(qubit_count, 0)
        sim  = similar(orig, shots=100)
        @test sim.shots  == 100
        @test orig.shots == 0
        sim.state_vector = 1/√2^qubit_count .* ones(ComplexF64, 2^qubit_count)
        sim2 = copy(sim)
        @test sim2.shots == 100
        @test sim2.state_vector  == 1/√2^qubit_count .* ones(ComplexF64, 2^qubit_count)
        sim.state_vector = 1/√2^qubit_count .* [exp(im*rand()) for ix in 1:2^qubit_count]
        copyto!(sim2, sim)
        @test sim2.state_vector ≈ sim.state_vector
    end
    @testset "supported operations and results" begin
        qubit_count = 10
        sim = StateVectorSimulator(qubit_count, 0)
        @test BraketSimulator.supported_operations(sim) == [
                "U",
                "GPhase",
                "ccnot",
                "cnot",
                "cphaseshift",
                "cphaseshift00",
                "cphaseshift01",
                "cphaseshift10",
                "cswap",
                "cv",
                "cy",
                "cz",
                "ecr",
                "gpi",
                "gpi2",
                "h",
                "i",
                "iswap",
                "ms",
                "pswap",
                "phaseshift",
                "prx",
                "rx",
                "ry",
                "rz",
                "s",
                "si",
                "swap",
                "t",
                "ti",
                "unitary",
                "v",
                "vi",
                "x",
                "xx",
                "xy",
                "y",
                "yy",
                "z",
                "zz",
            ]
    end
end
