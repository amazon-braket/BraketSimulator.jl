using Test, BraketSimulator, DataStructures, LinearAlgebra, BraketSimulator.Combinatorics

LARGE_TESTS = get(ENV, "BRAKET_SIM_LARGE_TESTS", "false") == "true"

@testset "State vector simulator" begin
    @testset "Simple circuits nq: $qubit_count $instructions" for (
        instructions,
        qubit_count,
        state_vector,
        probability_amplitudes,
    ) in [
        ([BraketSimulator.Instruction(BraketSimulator.H(), [0])], 1, [0.70710678, 0.70710678], [0.5, 0.5]),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0])], 1, [0, 1], [0, 1]),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0])], 2, [0, 0, 1, 0], [0, 0, 1, 0]),
        ([BraketSimulator.Instruction(BraketSimulator.Y(), [0])], 1, [0, im], [0, 1]),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.X(), [1])], 2, [0, 0, 0, 1], [0, 0, 0, 1]),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Z(), [0])], 1, [0, -1], [0, 1]),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1])],
            2,
            [0, 0, 0, 1],
            [0, 0, 0, 1],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.CY(), [0, 1])],
            2,
            [0, 0, 0, im],
            [0, 0, 0, 1],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.CZ(), [0, 1])],
            2,
            [0, 0, 1, 0],
            [0, 0, 1, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Swap(), [0, 1])],
            2,
            [0, 1, 0, 0],
            [0, 1, 0, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.ISwap(), [0, 1])],
            2,
            [0, im, 0, 0],
            [0, 1, 0, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Swap(), [0, 2])],
            3,
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
        ),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.S(), [0])], 1, [0, im], [0, 1]),
        ([BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Si(), [0])], 1, [0, -im], [0, 1]),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.T(), [0])],
            1,
            [0, 0.70710678 + 0.70710678 * im],
            [0, 1],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Ti(), [0])],
            1,
            [0, 0.70710678 - 0.70710678 * im],
            [0, 1],
        ),
        ([BraketSimulator.Instruction(BraketSimulator.V(), [0])], 1, [0.5 + 0.5 * im, 0.5 - 0.5 * im], [0.5, 0.5]),
        ([BraketSimulator.Instruction(BraketSimulator.Vi(), [0])], 1, [0.5 - 0.5 * im, 0.5 + 0.5 * im], [0.5, 0.5]),
        ([BraketSimulator.Instruction(BraketSimulator.I(), [0])], 1, [1, 0], [1, 0]),
        ([BraketSimulator.Instruction(BraketSimulator.Unitary([0 1; 1 0]), [0])], 1, [0, 1], [0, 1]),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.PhaseShift(0.15), [0])],
            1,
            [0, 0.98877108 + 0.14943813 * im],
            [0, 1],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.CPhaseShift(0.15), [0, 1]),
            ],
            2,
            [0, 0, 0, 0.98877108 + 0.14943813 * im],
            [0, 0, 0, 1],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.CPhaseShift00(0.15), [0, 1])],
            2,
            [0.98877108 + 0.14943813 * im, 0, 0, 0],
            [1, 0, 0, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [1]), BraketSimulator.Instruction(BraketSimulator.CPhaseShift01(0.15), [0, 1])],
            2,
            [0, 0.98877108 + 0.14943813 * im, 0, 0],
            [0, 1, 0, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.CPhaseShift10(0.15), [0, 1])],
            2,
            [0, 0, 0.98877108 + 0.14943813 * im, 0],
            [0, 0, 1, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.Rx(0.15), [0])],
            1,
            [0.99718882, -0.07492971 * im],
            [0.99438554, 0.00561446],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.Ry(0.15), [0])],
            1,
            [-0.07492971, 0.99718882],
            [0.00561446, 0.99438554],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.Rz(0.15), [0])],
            1,
            [0.70511898 - 0.0529833 * im, 0.70511898 + 0.0529833 * im],
            [0.5, 0.5],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.GPi(0.15), [0])],
            1,
            [0.69916673 - 0.10566872im, 0.69916673 + 0.10566872im],
            [0.5, 0.5],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.GPi2(0.15), [0])],
            1,
            [0.42528093 - 0.49438554im, 0.57471907 - 0.49438554im],
            [0.42528093, 0.57471907],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.MS(π/2, -π/4, 0.3), [0, 1])],
            2,
            [0.98877108, 0, 0, 0.10566872 - 0.10566872im],
            [0.97766824, 0, 0, 0.02233176],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1]), BraketSimulator.Instruction(BraketSimulator.ECR(), [0, 1])],
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
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.U(0.15, 0.4, 0.7), [0])],
            1,
            [0.66459511 - 0.03413278im, 0.36864009 + 0.64903989im],
            [0.44285171, 0.55714829],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.MultiQubitPhaseShift{1}(0.15), [0])],
            1,
            [0.69916673 + 0.10566872im, 0.69916673 + 0.10566872im],
            [0.5, 0.5],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.PRx(0.15, 0.4), [0])],
            1,
            [0.6844863 - 0.04880085im, 0.72575165 - 0.04880085im],
            [0.47090303, 0.52909697],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.PSwap(0.15), [0, 1])],
            2,
            [0, 0.98877108 + 0.14943813 * im, 0, 0],
            [0, 1, 0, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.X(), [0]), BraketSimulator.Instruction(BraketSimulator.XY(0.15), [0, 1])],
            2,
            [0, 0.07492971 * im, 0.99718882, 0],
            [0, 0.00561446, 0.99438554, 0],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.XX(0.3), [0, 1])],
            2,
            [0.98877108, 0, 0, -0.14943813 * im],
            [0.97766824, 0, 0, 0.02233176],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.YY(0.3), [0, 1])],
            2,
            [0.98877108, 0, 0, 0.14943813 * im],
            [0.97766824, 0, 0, 0.02233176],
        ),
        (
            [BraketSimulator.Instruction(BraketSimulator.ZZ(0.15), [0, 1])],
            2,
            [0.99718882 - 0.07492971 * im, 0, 0, 0],
            [1, 0, 0, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 1, 2]),
            ],
            3,
            [0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 1],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [0, 1, 2]),
            ],
            3,
            [0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [2]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [0, 2, 1]),
            ],
            3,
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.X(), [2]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [1, 2, 0]),
            ],
            3,
            [0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 1, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.X(), [2]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [2, 1, 0]),
            ],
            3,
            [0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [2]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [2, 0, 1]),
            ],
            3,
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
        ),
        (
            [
                BraketSimulator.Instruction(BraketSimulator.X(), [0]),
                BraketSimulator.Instruction(BraketSimulator.X(), [1]),
                BraketSimulator.Instruction(BraketSimulator.CSwap(), [1, 0, 2]),
            ],
            3,
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
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
        ([(BraketSimulator.Observables.X(), [0])], [BraketSimulator.Instruction(BraketSimulator.H(), [0])], 1),
        ([(BraketSimulator.Observables.Z(), [0])], BraketSimulator.Instruction[], 1),
        ([(BraketSimulator.Observables.I(), [0])], BraketSimulator.Instruction[], 1),
        (
            [
                (BraketSimulator.Observables.X(), [0]),
                (BraketSimulator.Observables.Z(), [3]),
                (BraketSimulator.Observables.H(), [2]),
            ],
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.Ry(-π / 4), [2])],
            5,
        ),
        (
            [(
                BraketSimulator.Observables.TensorProduct([
                    BraketSimulator.Observables.X(),
                    BraketSimulator.Observables.Z(),
                    BraketSimulator.Observables.H(),
                    BraketSimulator.Observables.I(),
                ]),
                (0, 3, 2, 1),
            )],
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.Ry(-π / 4), [2])],
            5,
        ),
        (
            [(BraketSimulator.Observables.X(), [0, 1])],
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.H(), [1])],
            2,
        ),
        ([(BraketSimulator.Observables.Z(), [0, 1])], BraketSimulator.Instruction[], 2),
        ([(BraketSimulator.Observables.I(), [0, 1])], BraketSimulator.Instruction[], 2),
        (
            [(
                BraketSimulator.Observables.TensorProduct([
                    BraketSimulator.Observables.I(),
                    BraketSimulator.Observables.Z(),
                ]),
                (2, 0),
            )],
            BraketSimulator.Instruction[],
            3,
        ),
        (
            [(
                BraketSimulator.Observables.TensorProduct([
                    BraketSimulator.Observables.X(),
                    BraketSimulator.Observables.Z(),
                ]),
                (2, 0),
            )],
            [BraketSimulator.Instruction(BraketSimulator.H(), [2])],
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
            [(BraketSimulator.Observables.X(), [0])],
        )
        @test_throws ErrorException BraketSimulator.apply_observables!(
            simulation,
            [(BraketSimulator.Observables.X(), [0])],
        )
    end
    @testset "state_with_observables fails before any observables are applied" begin
        simulation = StateVectorSimulator(4, 0)
        @test_throws ErrorException BraketSimulator.state_with_observables(simulation)
    end
    @testset "QFT simulation" begin
        function qft_circuit_operations(qubit_count::Int)
            qft_ops = BraketSimulator.Instruction[]
            for target_qubit = 0:qubit_count-1
                angle = π / 2
                push!(qft_ops, BraketSimulator.Instruction(BraketSimulator.H(), [target_qubit]))
                for control_qubit = target_qubit+1:qubit_count-1
                    push!(
                        qft_ops,
                        BraketSimulator.Instruction(BraketSimulator.CPhaseShift(angle), [control_qubit, target_qubit]),
                    )
                    angle /= 2
                end
            end
            return qft_ops
        end
        max_qc = 20
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
            [BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1])],
        )
        samples = counter(BraketSimulator.samples(simulation))

        @test BraketSimulator.qubit_count(simulation) == 2
        @test collect(keys(samples)) == [0, 3]
        @test 0.4 < samples[0] / (samples[0] + samples[3]) < 0.6
        @test 0.4 < samples[3] / (samples[0] + samples[3]) < 0.6
        @test samples[0] + samples[3] == 10000
    end
    @testset "batch" begin
        function make_ghz(num_qubits)
            ghz = BraketSimulator.Program(BraketSimulator.braketSchemaHeader("braket.ir.jaqcd.program", "1"), BraketSimulator.Instruction[BraketSimulator.Instruction(BraketSimulator.H(), 0)], BraketSimulator.AbstractProgramResult[], BraketSimulator.Instruction[])
            for ii in 1:num_qubits-1
                push!(ghz.instructions, BraketSimulator.Instruction(BraketSimulator.CNot(), [0, ii]))
            end
            return ghz
        end
        num_qubits = 5
        @testset for n_circuits in (1, 100)
            shots   = 1000
            jl_ghz  = [make_ghz(num_qubits) for ix in 1:n_circuits]
            jl_sim  = StateVectorSimulator(num_qubits, 0);
            results = BraketSimulator.simulate(jl_sim, jl_ghz, shots)
            for (r_ix, r) in enumerate(results)
                @test length(r.measurements) == shots
                @test 400 < count(m->m == fill(0, num_qubits), r.measurements) < 600
                @test 400 < count(m->m == fill(1, num_qubits), r.measurements) < 600
            end
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
    @testset "Measure with no gates" begin
        qasm = """
        bit[4] b;
        qubit[4] q;
        b[0] = measure q[0];
        b[1] = measure q[1];
        b[2] = measure q[2];
        b[3] = measure q[3];
        """
        simulator = StateVectorSimulator(0, 0)
        oq3_program = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), qasm, nothing)
        res = BraketSimulator.simulate(simulator, oq3_program, 1000)
        @test all(m -> m == zeros(Int, 4), res.measurements)
        @test length(res.measurements) == 1000
        @test res.measuredQubits == collect(0:3)
    end
    @testset "Measure with qubits not used" begin
        qasm = """
        bit[4] b;
        qubit[4] q;
        h q[0];
        cnot q[0], q[1];
        b = measure q;
        """ 
        simulator = StateVectorSimulator(0, 0)
        oq3_program = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), qasm, nothing)
        res = BraketSimulator.simulate(simulator, oq3_program, 1000)
        @test res.measuredQubits == collect(0:3)
        @test 400 < sum(m[1] for m in res.measurements) < 600
        @test 400 < sum(m[2] for m in res.measurements) < 600
        @test sum(m[3] for m in res.measurements) == 0
        @test sum(m[4] for m in res.measurements) == 0
        @test all(m->length(m) == 4, res.measurements)
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
        new_sv_qubit_count = 18
        new_sv_max_shots = 1_000_000
        new_sv_observables = ["x", "y", "z", "h", "i", "hermitian"]
        new_sv_props_dict = Dict(
            :braketSchemaHeader => Dict(
                :name => "braket.device_schema.simulators.gate_model_simulator_device_capabilities",
                :version => "1",
            ),
            :service => Dict(
                :executionWindows => [
                    Dict(
                        :executionDay => "Everyday",
                        :windowStartHour => "00:00",
                        :windowEndHour => "23:59:59",
                    ),
                ],
                :shotsRange => [0, new_sv_max_shots],
            ),
            :action => Dict(
                "braket.ir.jaqcd.program" => Dict(
                    :actionType => "braket.ir.jaqcd.program",
                    :version => ["1"],
                    :supportedOperations => [
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
                        "h",
                        "i",
                        "iswap",
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
                    ],
                    :supportedResultTypes => [
                        Dict(
                            :name => "Sample",
                            :observables => new_sv_observables,
                            :minShots => 1,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(
                            :name => "Expectation",
                            :observables => new_sv_observables,
                            :minShots => 0,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(
                            :name => "Variance",
                            :observables => new_sv_observables,
                            :minShots => 0,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(:name => "Probability", :minShots => 0, :maxShots => new_sv_max_shots),
                        Dict(:name => "StateVector", :minShots => 0, :maxShots => 0),
                        Dict(:name => "DensityMatrix", :minShots => 0, :maxShots => 0),
                        Dict(:name => "Amplitude", :minShots => 0, :maxShots => 0),
                    ],
                ),
                "braket.ir.openqasm.program" => Dict(
                    :actionType => "braket.ir.openqasm.program",
                    :version => ["1"],
                    :supportedOperations => [
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
                    ],
                    :supportedModifiers => [
                        Dict(:name => "ctrl"),
                        Dict(:name => "negctrl"),
                        Dict(:name => "pow", :exponent_types => ["int", "float"]),
                        Dict(:name => "inv"),
                    ],
                    :supportedPragmas => [
                        "braket_unitary_matrix",
                        "braket_result_type_state_vector",
                        "braket_result_type_density_matrix",
                        "braket_result_type_sample",
                        "braket_result_type_expectation",
                        "braket_result_type_variance",
                        "braket_result_type_probability",
                        "braket_result_type_amplitude",
                    ],
                    :forbiddenPragmas => [
                        "braket_noise_amplitude_damping",
                        "braket_noise_bit_flip",
                        "braket_noise_depolarizing",
                        "braket_noise_kraus",
                        "braket_noise_pauli_channel",
                        "braket_noise_generalized_amplitude_damping",
                        "braket_noise_phase_flip",
                        "braket_noise_phase_damping",
                        "braket_noise_two_qubit_dephasing",
                        "braket_noise_two_qubit_depolarizing",
                    ],
                    :supportedResultTypes => [
                        Dict(
                            :name => "Sample",
                            :observables => new_sv_observables,
                            :minShots => 1,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(
                            :name => "Expectation",
                            :observables => new_sv_observables,
                            :minShots => 0,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(
                            :name => "Variance",
                            :observables => new_sv_observables,
                            :minShots => 0,
                            :maxShots => new_sv_max_shots,
                        ),
                        Dict(:name => "Probability", :minShots => 0, :maxShots => new_sv_max_shots),
                        Dict(:name => "StateVector", :minShots => 0, :maxShots => 0),
                        Dict(:name => "DensityMatrix", :minShots => 0, :maxShots => 0),
                        Dict(:name => "Amplitude", :minShots => 0, :maxShots => 0),
                    ],
                    :supportPhysicalQubits => false,
                    :supportsPartialVerbatimBox => false,
                    :requiresContiguousQubitIndices => true,
                    :requiresAllQubitsMeasurement => true,
                    :supportsUnassignedMeasurements => true,
                    :disabledQubitRewiringSupported => false,
                ),
            ),
            :paradigm => Dict(:qubitCount => new_sv_qubit_count),
            :deviceParameters =>
                Dict(:paradigmParameters => Dict(:qubitCount => new_sv_qubit_count)),
        )

        new_sv_props = BraketSimulator.StructTypes.constructfrom(BraketSimulator.GateModelSimulatorDeviceCapabilities, new_sv_props_dict)
        @test new_sv_props.paradigm.qubitCount == new_sv_qubit_count
        @test BraketSimulator.supported_result_types(sim) == BraketSimulator.supported_result_types(sim, Val(:OpenQASM))
    end
    @testset "mmaping large results" begin
        oq3_source = """OPENQASM 3.0;\nqubit[18] q;\nh q;\n#pragma braket result state_vector\n"""
        oq3_result  = BraketSimulator.simulate("braket_sv_v2", oq3_source, "{}", 0)
        result      = JSON3.read(oq3_result[1], BraketSimulator.GateModelTaskResult)
        @test isempty(result.resultTypes[1].value)
        @test !isnothing(oq3_result[2]) && !isempty(oq3_result[2])
        @test !isnothing(oq3_result[3]) && !isempty(oq3_result[3])
    end
    @testset "partial trace $nq" for nq in 3:6
        ψ       = normalize(rand(ComplexF64, 2^nq))
        full_ρ  = kron(ψ, adjoint(ψ))
        @testset "output qubits $q" for q in combinations(0:nq-1)
            ρ       = BraketSimulator.partial_trace(ψ, q)
            full_pt = BraketSimulator.partial_trace(full_ρ, q)
            @test full_pt ≈ ρ
        end
    end
    @testset "inputs handling" begin
        sv_adder_qasm = """
            OPENQASM 3;

            input uint[4] a_in;
            input uint[4] b_in;

            gate majority a, b, c {
                cnot c, b;
                cnot c, a;
                ccnot a, b, c;
            }

            gate unmaj a, b, c {
                ccnot a, b, c;
                cnot c, a;
                cnot a, b;
            }

            qubit cin;
            qubit[4] a;
            qubit[4] b;
            qubit cout;

            // set input states
            for int[8] i in [0: 3] {
              if(bool(a_in[i])) x a[i];
              if(bool(b_in[i])) x b[i];
            }

            // add a to b, storing result in b
            majority cin, b[3], a[3];
            for int[8] i in [3: -1: 1] { majority a[i], b[i - 1], a[i - 1]; }
            cnot a[0], cout;
            for int[8] i in [1: 3] { unmaj a[i], b[i - 1], a[i - 1]; }
            unmaj cin, b[3], a[3];

            // todo: subtle bug when trying to get a result type for both at once
            #pragma braket result probability cout, b
            #pragma braket result probability cout
            #pragma braket result probability b
            """
        simulator = StateVectorSimulator(6, 0)
        oq3_program = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), sv_adder_qasm, Dict("a_in"=>3, "b_in"=>7))
        # should NOT throw a missing input error
        BraketSimulator.simulate(simulator, oq3_program, 0; inputs=Dict{String, Float64}())
    end
end
