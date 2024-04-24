using Test, PythonCall, BraketSimulator, Braket

@testset "Python integration" begin
    braket_ixs       = pyimport("braket.ir.jaqcd.instructions")
    braket_rts       = pyimport("braket.ir.jaqcd.results")
    braket_sim_gates = pyimport("braket.default_simulator.gate_operations")
    braket_sim_noise = pyimport("braket.default_simulator.noise_operations")
    braket_sim_rts   = pyimport("braket.default_simulator.result_types")
    braket_sim_obs   = pyimport("braket.default_simulator.observables")
    np               = pyimport("numpy")
    jl_mat           = ComplexF64[0.0 1.0; 1.0 0.0]
    py_mat           = pylist([pylist([pylist([0.0; 0.0]); pylist([1.0; 0.0])]); pylist([pylist([1.0; 0.0]); pylist([0.0; 0.0])])])
    py_sim_mat       = np.array(pylist([pylist([0.0, 1.0]); pylist([1.0, 0.0])]))
    @testset "Gates" begin
        @testset "1q static gate $jl_gate" for (jl_gate, py_gate, py_sim_gate) in zip(
            [Braket.I(), H(), X(), Y(), Z(), V(), Vi(), T(), Ti(), S(), Si(), Unitary(jl_mat)],
            [
                braket_ixs.I(target=0),
                braket_ixs.H(target=0),
                braket_ixs.X(target=0),
                braket_ixs.Y(target=0),
                braket_ixs.Z(target=0),
                braket_ixs.V(target=0),
                braket_ixs.Vi(target=0),
                braket_ixs.T(target=0),
                braket_ixs.Ti(target=0),
                braket_ixs.S(target=0),
                braket_ixs.Si(target=0),
                braket_ixs.Unitary(targets=pylist([0]), matrix=py_mat),
            ],
            [
                braket_sim_gates.Identity(targets=pylist([0])),
                braket_sim_gates.Hadamard(targets=pylist([0])),
                braket_sim_gates.PauliX(targets=pylist([0])),
                braket_sim_gates.PauliY(targets=pylist([0])),
                braket_sim_gates.PauliZ(targets=pylist([0])),
                braket_sim_gates.V(targets=pylist([0])),
                braket_sim_gates.Vi(targets=pylist([0])),
                braket_sim_gates.T(targets=pylist([0])),
                braket_sim_gates.Ti(targets=pylist([0])),
                braket_sim_gates.S(targets=pylist([0])),
                braket_sim_gates.Si(targets=pylist([0])),
                braket_sim_gates.Unitary(targets=pylist([0]), matrix=py_sim_mat),
            ],
        )
            @test pyconvert(Bool, Py(jl_gate, 0) == py_gate)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_gate, [0])) == py_gate)
            @test pyconvert(typeof(jl_gate), py_gate) == jl_gate
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        angle    = π / 3.5
        angle2   = π / 5.2
        angle3   = π / 0.6
        @testset "1q 1-parameter gate $jl_gate" for (jl_gate, py_gate, py_sim_gate) in zip(
            [Rx(angle), Ry(angle), Rz(angle), PhaseShift(angle)],
            [
                braket_ixs.Rx(target=0, angle=angle),
                braket_ixs.Ry(target=0, angle=angle),
                braket_ixs.Rz(target=0, angle=angle),
                braket_ixs.PhaseShift(target=0, angle=angle),
            ],
            [
                braket_sim_gates.RotX(targets=pylist([0]), angle=angle),
                braket_sim_gates.RotY(targets=pylist([0]), angle=angle),
                braket_sim_gates.RotZ(targets=pylist([0]), angle=angle),
                braket_sim_gates.PhaseShift(targets=pylist([0]), angle=angle),
            ],
        )
            @test pyconvert(Bool, Py(jl_gate, 0) == py_gate)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_gate, [0])) == py_gate)
            @test pyconvert(typeof(jl_gate), py_gate) == jl_gate
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        @testset "1q 1-parameter gate $jl_gate" for (jl_gate, py_sim_gate) in zip(
            [
                GPi(angle),
                GPi2(angle),
                MultiQubitPhaseShift{1}(angle),
            ],
            [
                braket_sim_gates.GPi(targets=pylist([0]), angle=angle),
                braket_sim_gates.GPi2(targets=pylist([0]), angle=angle),
                braket_sim_gates.GPhase(targets=pylist([0]), angle=angle),
            ],
        )
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        @testset "2-parameter gate $jl_gate" for (jl_gate, py_sim_gate, targets) in zip(
            [
                PRx(angle, angle2),
            ],
            [
                braket_sim_gates.PRx(targets=pylist([0]), angle_1=angle, angle_2=angle2),
            ],
            [
                [0],
            ],
        )
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, targets)
        end
        @testset "3-parameter gate $jl_gate" for (jl_gate, py_sim_gate, targets) in zip(
            [
                MS(angle, angle2, angle3),
                U(angle, angle2, angle3),
            ],
            [
                braket_sim_gates.MS(targets=pylist([0, 1]), angle_1=angle, angle_2=angle2, angle_3=angle3),
                braket_sim_gates.U(targets=pylist([0]), theta=angle, phi=angle2, lambda_=angle3, ctrl_modifiers=pylist()),
            ],
            [
                [0, 1],
                [0],
            ],
        )
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, targets)
        end
        @testset "2q static gate $jl_gate" for (jl_gate, py_gate, py_sim_gate) in zip(
            [CNot(), CY(), CZ(), CV(), Swap(), ISwap(), ECR()],
            [
                braket_ixs.CNot(control=0, target=1),
                braket_ixs.CY(control=0, target=1),
                braket_ixs.CZ(control=0, target=1),
                braket_ixs.CV(control=0, target=1),
                braket_ixs.Swap(targets=pylist([0, 1])),
                braket_ixs.ISwap(targets=pylist([0, 1])),
                braket_ixs.ECR(targets=pylist([0, 1])),
            ],
            [
                braket_sim_gates.CX(targets=pylist([0, 1])),
                braket_sim_gates.CY(targets=pylist([0, 1])),
                braket_sim_gates.CZ(targets=pylist([0, 1])),
                braket_sim_gates.CV(targets=pylist([0, 1])),
                braket_sim_gates.Swap(targets=pylist([0, 1])),
                braket_sim_gates.ISwap(targets=pylist([0, 1])),
                braket_sim_gates.ECR(targets=pylist([0, 1])),
            ],
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1]) == py_gate)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_gate, [0, 1])) == py_gate)
            @test pyconvert(typeof(jl_gate), py_gate) == jl_gate
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset "2q 1-parameter gate $jl_gate" for (jl_gate, py_gate, py_sim_gate) in zip(
            [
                XX(angle),
                XY(angle),
                YY(angle),
                ZZ(angle),
                CPhaseShift(angle),
                CPhaseShift00(angle),
                CPhaseShift10(angle),
                CPhaseShift01(angle),
                PSwap(angle),
            ],
            [
                braket_ixs.XX(targets=pylist([0, 1]), angle=angle),
                braket_ixs.XY(targets=pylist([0, 1]), angle=angle),
                braket_ixs.YY(targets=pylist([0, 1]), angle=angle),
                braket_ixs.ZZ(targets=pylist([0, 1]), angle=angle),
                braket_ixs.CPhaseShift(control=0, target=1, angle=angle),
                braket_ixs.CPhaseShift00(control=0, target=1, angle=angle),
                braket_ixs.CPhaseShift10(control=0, target=1, angle=angle),
                braket_ixs.CPhaseShift01(control=0, target=1, angle=angle),
                braket_ixs.PSwap(targets=pylist([0, 1]), angle=angle),
            ],
            [
                braket_sim_gates.XX(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.XY(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.YY(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.ZZ(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.CPhaseShift(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.CPhaseShift00(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.CPhaseShift10(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.CPhaseShift01(targets=pylist([0, 1]), angle=angle),
                braket_sim_gates.PSwap(targets=pylist([0, 1]), angle=angle),
            ],
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1]) == py_gate)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_gate, [0, 1])) == py_gate)
            @test pyconvert(typeof(jl_gate), py_gate) == jl_gate
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset "3q static gate $jl_gate" for (jl_gate, py_gate, py_sim_gate) in zip(
            [CCNot(), CSwap()],
            [braket_ixs.CCNot(controls=pylist([0, 1]), target=2),
             braket_ixs.CSwap(control=0, targets=pylist([1, 2]))],
            [braket_sim_gates.CCNot(targets=pylist([0, 1, 2])),
             braket_sim_gates.CSwap(targets=pylist([0, 1, 2]))],
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1, 2]) == py_gate)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_gate, [0, 1, 2])) == py_gate)
            @test pyconvert(typeof(jl_gate), py_gate) == jl_gate
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1, 2])
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, [0, 1, 2])
        end
    end
    prob = 0.15
    gamma = 0.25
    probx = 0.1
    proby = 0.2
    probz = 0.3
    @testset "Noises" begin
        @testset "1q noise $jl_noise" for (jl_noise, py_noise, py_sim_noise) in zip(
                [
                    BitFlip(prob),
                    PhaseFlip(prob),
                    Depolarizing(prob),
                    AmplitudeDamping(gamma),
                    GeneralizedAmplitudeDamping(prob, gamma),
                    PhaseDamping(gamma),
                    PauliChannel(probx, proby, probz),
                    Kraus([jl_mat]),
                ],
                [
                    braket_ixs.BitFlip(target=0, probability=prob),
                    braket_ixs.PhaseFlip(target=0, probability=prob),
                    braket_ixs.Depolarizing(target=0, probability=prob),
                    braket_ixs.AmplitudeDamping(target=0, gamma=gamma),
                    braket_ixs.GeneralizedAmplitudeDamping(target=0, probability=prob, gamma=gamma),
                    braket_ixs.PhaseDamping(target=0, gamma=gamma),
                    braket_ixs.PauliChannel(target=0, probX=probx, probY=proby, probZ=probz),
                    braket_ixs.Kraus(targets=pylist([0]), matrices=pylist([py_mat])),
                ],
                [
                    braket_sim_noise.BitFlip(targets=pylist([0]), probability=prob),
                    braket_sim_noise.PhaseFlip(targets=pylist([0]), probability=prob),
                    braket_sim_noise.Depolarizing(targets=pylist([0]), probability=prob),
                    braket_sim_noise.AmplitudeDamping(targets=pylist([0]), gamma=gamma),
                    braket_sim_noise.GeneralizedAmplitudeDamping(targets=pylist([0]), probability=prob, gamma=gamma),
                    braket_sim_noise.PhaseDamping(targets=pylist([0]), gamma=gamma),
                    braket_sim_noise.PauliChannel(targets=pylist([0]), probX=probx, probY=proby, probZ=probz),
                    braket_sim_noise.Kraus(targets=pylist([0]), matrices=pylist([py_sim_mat])),
                ],
            )
            @test pyconvert(Bool, Py(jl_noise, 0) == py_noise)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_noise, [0])) == py_noise)
            @test pyconvert(typeof(jl_noise), py_noise) == jl_noise
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, 0)
            @test pyconvert(Braket.Instruction, py_sim_noise) == Braket.Instruction(jl_noise, 0)
        end
        @testset "2q noise $jl_noise" for (jl_noise, py_noise, py_sim_noise) in zip(
            [TwoQubitDepolarizing(prob), TwoQubitDephasing(prob)],
            [
             braket_ixs.TwoQubitDepolarizing(targets=pylist([0, 1]), probability=prob),
             braket_ixs.TwoQubitDephasing(targets=pylist([0, 1]), probability=prob),
            ],
            [
             braket_sim_noise.TwoQubitDepolarizing(targets=pylist([0, 1]), probability=prob),
             braket_sim_noise.TwoQubitDephasing(targets=pylist([0, 1]), probability=prob),
            ],
        )
            @test pyconvert(Bool, Py(jl_noise, [0, 1]) == py_noise)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_noise, [0, 1])) == py_noise)
            @test pyconvert(typeof(jl_noise), py_noise) == jl_noise
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_noise) == Braket.Instruction(jl_noise, [0, 1])
        end
    end
    @testset "Compiler directives" begin
        @testset for (jl_cd, py_cd) in zip([StartVerbatimBox(), EndVerbatimBox()],
                                           [braket_ixs.StartVerbatimBox(), braket_ixs.EndVerbatimBox()])
            @test pyconvert(Bool, Py(jl_cd) == py_cd)
            @test pyconvert(Bool, Py(Braket.Instruction(jl_cd)) == py_cd)
            @test pyconvert(typeof(jl_cd), py_cd)      == jl_cd
            @test pyconvert(Braket.Instruction, py_cd) == Braket.Instruction(jl_cd)
        end
    end
    @testset "Result types" begin
        jl_tp_1 = ["h", "x"]
        py_tp_1 = pylist([pystr("h"), pystr("x")])
        jl_herm = Braket.complex_matrix_to_ir(jl_mat)
        py_herm = py_mat
        jl_tp_2 = Union{String, Vector{Vector{Vector{Float64}}}}["h", jl_herm]
        py_tp_2 = pylist([pystr("h"), py_herm])
        @testset for (jl_rt, py_rt, jl_braket_rt) in zip(
                [
                    Braket.IR.StateVector("statevector"),
                    Braket.IR.Amplitude(["00", "11"], "amplitude"),
                    Braket.IR.Probability(nothing, "probability"),
                    Braket.IR.Probability([0], "probability"),
                    Braket.IR.DensityMatrix(nothing, "densitymatrix"),
                    Braket.IR.DensityMatrix([0], "densitymatrix"),
                    Braket.IR.Sample("x", nothing, "sample"),
                    Braket.IR.Sample("x", [0], "sample"),
                    Braket.IR.Sample([jl_herm], nothing, "sample"),
                    Braket.IR.Sample([jl_herm], [0], "sample"),
                    Braket.IR.Sample(jl_tp_1, [0, 1], "sample"),
                    Braket.IR.Sample(jl_tp_2, [0, 1], "sample"),
                    Braket.IR.Expectation("x", nothing, "expectation"),
                    Braket.IR.Expectation("x", [0], "expectation"),
                    Braket.IR.Expectation([jl_herm], nothing, "expectation"),
                    Braket.IR.Expectation([jl_herm], [0], "expectation"),
                    Braket.IR.Expectation(jl_tp_1, [0, 1], "expectation"),
                    Braket.IR.Expectation(jl_tp_2, [0, 1], "expectation"),
                    Braket.IR.Variance("x", nothing, "variance"),
                    Braket.IR.Variance("x", [0], "variance"),
                    Braket.IR.Variance([jl_herm], nothing, "variance"),
                    Braket.IR.Variance([jl_herm], [0], "variance"),
                    Braket.IR.Variance(jl_tp_1, [0, 1], "variance"),
                    Braket.IR.Variance(jl_tp_2, [0, 1], "variance"),
                ],
                [
                    braket_rts.StateVector(),
                    braket_rts.Amplitude(states=pylist([pystr("00"), pystr("11")])),
                    braket_rts.Probability(),
                    braket_rts.Probability(targets=pylist([0])),
                    braket_rts.DensityMatrix(),
                    braket_rts.DensityMatrix(targets=pylist([0])),
                    braket_rts.Sample(observable=pylist([pystr("x")])),
                    braket_rts.Sample(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Sample(observable=pylist([py_herm])),
                    braket_rts.Sample(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Sample(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Sample(observable=py_tp_2, targets=pylist([0, 1])),
                    braket_rts.Expectation(observable=pylist([pystr("x")])),
                    braket_rts.Expectation(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Expectation(observable=pylist([py_herm])),
                    braket_rts.Expectation(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Expectation(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Expectation(observable=py_tp_2, targets=pylist([0, 1])),
                    braket_rts.Variance(observable=pylist([pystr("x")])),
                    braket_rts.Variance(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Variance(observable=pylist([py_herm])),
                    braket_rts.Variance(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Variance(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Variance(observable=py_tp_2, targets=pylist([0, 1])),
                ],
                [
                    Braket.StateVector(),
                    Braket.Amplitude(["00", "11"]),
                    Braket.Probability(),
                    Braket.Probability(0),
                    Braket.DensityMatrix(),
                    Braket.DensityMatrix(0),
                    Braket.Sample(Braket.Observables.X(), Int[]),
                    Braket.Sample(Braket.Observables.X(), [0]),
                    Braket.Sample(Braket.Observables.HermitianObservable(jl_mat), Int[]),
                    Braket.Sample(Braket.Observables.HermitianObservable(jl_mat), [0]),
                    Braket.Sample(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.X()]), [0, 1]),
                    Braket.Sample(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.HermitianObservable(jl_mat)]), [0, 1]),
                    Braket.Expectation(Braket.Observables.X(), Int[]),
                    Braket.Expectation(Braket.Observables.X(), [0]),
                    Braket.Expectation(Braket.Observables.HermitianObservable(jl_mat), Int[]),
                    Braket.Expectation(Braket.Observables.HermitianObservable(jl_mat), [0]),
                    Braket.Expectation(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.X()]), [0, 1]),
                    Braket.Expectation(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.HermitianObservable(jl_mat)]), [0, 1]),
                    Braket.Variance(Braket.Observables.X(), Int[]),
                    Braket.Variance(Braket.Observables.X(), [0]),
                    Braket.Variance(Braket.Observables.HermitianObservable(jl_mat), Int[]),
                    Braket.Variance(Braket.Observables.HermitianObservable(jl_mat), [0]),
                    Braket.Variance(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.X()]), [0, 1]),
                    Braket.Variance(Braket.Observables.TensorProduct([Braket.Observables.H(), Braket.Observables.HermitianObservable(jl_mat)]), [0, 1]),
                ],
        )
            @test pyconvert(Bool, Py(jl_rt) == py_rt)
            @test pyconvert(Braket.AbstractProgramResult, py_rt) == jl_rt
            @test pyconvert(typeof(jl_braket_rt), py_rt) == jl_braket_rt
            @test pyconvert(Result, py_rt) == jl_braket_rt
        end
    end
    @testset "Result types and observables from default sim" begin
        jl_tp_1 = Braket.Observables.TensorProduct(["h", "x"])
        py_tp_1 = braket_sim_obs.TensorProduct(pylist([braket_sim_obs.Hadamard(targets=pylist([0])), braket_sim_obs.PauliX(targets=pylist([1]))]))
        jl_herm = Braket.Observables.HermitianObservable(jl_mat)
        py_herm = braket_sim_obs.Hermitian(targets=pylist([0]), matrix=py_sim_mat)
        jl_tp_2 = Braket.Observables.TensorProduct([Braket.Observables.H(), jl_herm])
        py_tp_2 = braket_sim_obs.TensorProduct(pylist([braket_sim_obs.Hadamard(targets=pylist([0])), braket_sim_obs.Hermitian(targets=pylist([1]), matrix=py_sim_mat)]))
        @testset for (jl_obs, py_obs) in zip(
                [
                    jl_tp_1,
                    jl_tp_2,
                    jl_herm,
                    Braket.Observables.H(),
                    Braket.Observables.X(),
                    Braket.Observables.Y(),
                    Braket.Observables.Z(),
                    Braket.Observables.I(),
                ],
                [
                    py_tp_1,
                    py_tp_2,
                    py_herm,
                    braket_sim_obs.Hadamard(targets=pylist([0])),
                    braket_sim_obs.PauliX(targets=pylist([0])),
                    braket_sim_obs.PauliY(targets=pylist([0])),
                    braket_sim_obs.PauliZ(targets=pylist([0])),
                    braket_sim_obs.Identity(targets=pylist([0])),
                ],
               )
            # do this to avoid overspecialized parameter
            T = jl_obs isa Braket.Observables.TensorProduct ? Braket.Observables.TensorProduct : typeof(jl_obs)
            @test pyconvert(T, py_obs) == jl_obs
            @test pyconvert(Braket.Observables.Observable, py_obs) == jl_obs
        end
        @testset for (jl_rt, py_rt) in zip(
                [
                    Braket.StateVector(),
                    Braket.Amplitude(["00", "11"]),
                    Braket.Probability(QubitSet()),
                    Braket.Probability([0]),
                    Braket.DensityMatrix(QubitSet()),
                    Braket.DensityMatrix([0]),
                    Braket.Expectation(Braket.Observables.X(), [0]),
                    Braket.Expectation(jl_herm, [0]),
                    Braket.Expectation(jl_tp_1, [0, 1]),
                    Braket.Expectation(jl_tp_2, [0, 1]),
                    Braket.Variance(Braket.Observables.X(), [0]),
                    Braket.Variance(jl_herm, [0]),
                    Braket.Variance(jl_tp_1, [0, 1]),
                    Braket.Variance(jl_tp_2, [0, 1]),
                ],
                [
                    braket_sim_rts.StateVector(),
                    braket_sim_rts.Amplitude(states=pylist([pystr("00"), pystr("11")])),
                    braket_sim_rts.Probability(),
                    braket_sim_rts.Probability(targets=pylist([0])),
                    braket_sim_rts.DensityMatrix(),
                    braket_sim_rts.DensityMatrix(targets=pylist([0])),
                    braket_sim_rts.Expectation(observable=braket_sim_obs.PauliX(targets=pylist([0]))),
                    braket_sim_rts.Expectation(observable=py_herm),
                    braket_sim_rts.Expectation(observable=py_tp_1),
                    braket_sim_rts.Expectation(observable=py_tp_2),
                    braket_sim_rts.Variance(observable=braket_sim_obs.PauliX(targets=pylist([0]))),
                    braket_sim_rts.Variance(observable=py_herm),
                    braket_sim_rts.Variance(observable=py_tp_1),
                    braket_sim_rts.Variance(observable=py_tp_2),
                ],
        )
            @test pyconvert(typeof(jl_rt), py_rt) == jl_rt
            @test pyconvert(Result, py_rt) == jl_rt
        end
        py_tp_1 = pylist([pystr("h"), pystr("x")])
        py_herm = py_mat
        py_tp_2 = pylist([pystr("h"), py_herm])
        @testset for (jl_rt, py_rt) in zip(
                [
                    Braket.StateVector(),
                    Braket.Amplitude(["00", "11"]),
                    Braket.Probability(),
                    Braket.Probability([0]),
                    Braket.DensityMatrix(),
                    Braket.DensityMatrix([0]),
                    Braket.Sample(Braket.Observables.X(), nothing),
                    Braket.Sample(Braket.Observables.X(), [0]),
                    Braket.Sample(jl_herm, nothing),
                    Braket.Sample(jl_herm, [0]),
                    Braket.Sample(jl_tp_1, [0, 1]),
                    Braket.Sample(jl_tp_2, [0, 1]),
                    Braket.Expectation(Braket.Observables.X(), nothing),
                    Braket.Expectation(Braket.Observables.X(), [0]),
                    Braket.Expectation(jl_herm, nothing),
                    Braket.Expectation(jl_herm, [0]),
                    Braket.Expectation(jl_tp_1, [0, 1]),
                    Braket.Expectation(jl_tp_2, [0, 1]),
                    Braket.Variance(Braket.Observables.X(), nothing),
                    Braket.Variance(Braket.Observables.X(), [0]),
                    Braket.Variance(jl_herm, nothing),
                    Braket.Variance(jl_herm, [0]),
                    Braket.Variance(jl_tp_1, [0, 1]),
                    Braket.Variance(jl_tp_2, [0, 1]),
                ],
                [
                    braket_rts.StateVector(),
                    braket_rts.Amplitude(states=pylist([pystr("00"), pystr("11")])),
                    braket_rts.Probability(),
                    braket_rts.Probability(targets=pylist([0])),
                    braket_rts.DensityMatrix(),
                    braket_rts.DensityMatrix(targets=pylist([0])),
                    braket_rts.Sample(observable=pylist([pystr("x")])),
                    braket_rts.Sample(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Sample(observable=pylist([py_herm])),
                    braket_rts.Sample(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Sample(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Sample(observable=py_tp_2, targets=pylist([0, 1])),
                    braket_rts.Expectation(observable=pylist([pystr("x")])),
                    braket_rts.Expectation(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Expectation(observable=pylist([py_herm])),
                    braket_rts.Expectation(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Expectation(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Expectation(observable=py_tp_2, targets=pylist([0, 1])),
                    braket_rts.Variance(observable=pylist([pystr("x")])),
                    braket_rts.Variance(observable=pylist([pystr("x")]), targets=pylist([0])),
                    braket_rts.Variance(observable=pylist([py_herm])),
                    braket_rts.Variance(observable=pylist([py_herm]), targets=pylist([0])),
                    braket_rts.Variance(observable=py_tp_1, targets=pylist([0, 1])),
                    braket_rts.Variance(observable=py_tp_2, targets=pylist([0, 1])),
                ],
        )
            @test pyconvert(typeof(jl_rt), py_rt) == jl_rt
            @test pyconvert(Result, py_rt) == jl_rt
        end
    end
    @testset "Programs" begin
        sv_adder = """
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
        oq3_program = Braket.OpenQasmProgram(Braket.braketSchemaHeader("braket.ir.openqasm.program", "1"), sv_adder, Dict("a_in"=>3, "b_in"=>7))
        n_qubits = 5
        function qft_circuit(qubit_count::Int)
            qft_circ = Circuit() 
            for target_qubit = 0:qubit_count-1
                angle = π/2
                qft_circ(H(), target_qubit)
                for control_qubit = target_qubit+1:qubit_count-1
                    qft_circ(CPhaseShift(angle), control_qubit, target_qubit)
                    angle /= 2
                end
            end
            qft_circ(Amplitude([repeat("0", qubit_count), repeat("1", qubit_count)]))
            qft_circ(Expectation(Braket.Observables.X(), 0))
            qft_circ(DensityMatrix())
            return qft_circ
        end
        jaqcd_program = ir(qft_circuit(n_qubits), Val(:JAQCD))
        @test pyconvert(Braket.OpenQasmProgram, Py(oq3_program)) == oq3_program
        @test pyconvert(Braket.Program, Py(jaqcd_program)) == jaqcd_program
        @testset "Full Python circuit execution" begin
            @testset "OpenQASM3" begin
                sv_simulator = StateVectorSimulator(n_qubits, 0)
                oq3_results  = simulate(sv_simulator, PyList{Any}([Py(oq3_program)]), 0)
                @test pyconvert(Bool, oq3_results.resultTypes[0].type == Py(Braket.IR.Probability([9, 5, 6, 7, 8], "probability")))
                # test a "batch"
                oq3_results  = simulate(sv_simulator, PyList{Any}([Py(oq3_program)]), 0; input = pylist([pydict(Dict("a_in"=>2, "b_in"=>5)), pydict(Dict("a_in"=>3, "b_in"=>2))]))
                @test pyconvert(Vector{Float64}, oq3_results.resultTypes[0].value) ≠ pyconvert(Vector{Float64}, oq3_results.resultTypes[1].value)
            end
            @testset "JAQCD" begin
                sv_simulator  = StateVectorSimulator(n_qubits, 0)
                jaqcd_results = simulate(sv_simulator, PyList{Any}([Py(jaqcd_program)]), n_qubits, 0)
                @test pyconvert(Bool, jaqcd_results.resultTypes[0].type == Py(Braket.IR.Amplitude([repeat("0", n_qubits), repeat("1", n_qubits)], "amplitude")))
                @test pyconvert(Bool, jaqcd_results.resultTypes[1].type == Py(Braket.IR.Expectation(["x"], [0], "expectation")))
                @test pyconvert(Bool, jaqcd_results.resultTypes[2].type == Py(Braket.IR.DensityMatrix(nothing, "density_matrix")))
            end
        end
    end
    @testset "Python circuit with measured qubits" begin
        qasm = """
        qubit[2] q; 
        bit[1] b;
        h q[0];
        cnot q[0], q[1];
        b[0] = measure q[0];
        """
        simulator    = StateVectorSimulator(2, 1000)
        oq3_program  = Braket.OpenQasmProgram(Braket.braketSchemaHeader("braket.ir.openqasm.program", "1"), qasm, nothing)
        result       = simulate(simulator, PyList{Any}([Py(oq3_program)]), 1000; measured_qubits=pylist([0]))
        measurements = [[pyconvert(Int, m) for m in measurement] for measurement in result.measurements]
        @test 400 < sum(m[1] for m in measurements) < 600
        @test all(length(m) == 1 for m in measurements)
        @test pyconvert(Bool, result.measuredQubits == pylist([0]))
    end
end
