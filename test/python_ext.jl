using Test, PythonCall, BraketSimulator, Braket

@testset "Python integration" begin
    braket_ixs = pyimport("braket.ir.jaqcd.instructions")
    braket_rts = pyimport("braket.ir.jaqcd.results")
    braket_sim_gates = pyimport("braket.default_simulator.gate_operations")
    braket_sim_noise = pyimport("braket.default_simulator.noise_operations")
    braket_sim_rts   = pyimport("braket.default_simulator.result_types")
    braket_sim_obs   = pyimport("braket.default_simulator.observables")
    np               = pyimport("numpy")
    jl_mat     = ComplexF64[0.0 1.0; 1.0 0.0]
    py_mat     = pylist([pylist([pylist([0.0; 0.0]); pylist([1.0; 0.0])]); pylist([pylist([1.0; 0.0]); pylist([0.0; 0.0])])])
    py_sim_mat = np.array(pylist([pylist([0.0, 1.0]); pylist([1.0, 0.0])]))
    @testset "Gates" begin
        @testset for (jl_gate, py_gate, py_sim_gate) in zip(
            [H(), X(), Y(), Z(), V(), Vi(), T(), Ti(), S(), Si(), Unitary(jl_mat)],
            [
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
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        angle    = π / 3.5
        angle2   = π / 5.2
        angle3   = π / 0.6
        @testset for (jl_gate, py_gate, py_sim_gate) in zip(
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
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        @testset for (jl_gate, py_sim_gate) in zip(
            [
                GPi(angle),
                GPi2(angle),
                #MultiQubitPhaseShift{1}(angle),
            ],
            [
                braket_sim_gates.GPi(targets=pylist([0]), angle=angle),
                braket_sim_gates.GPi2(targets=pylist([0]), angle=angle),
                #braket_sim_gates.GPhase(angle=angle),
            ],
        )
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, 0)
        end
        @testset for (jl_gate, py_sim_gate, targets) in zip(
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
        @testset for (jl_gate, py_gate, py_sim_gate) in zip(
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
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset for (jl_gate, py_gate, py_sim_gate) in zip(
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
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset for (jl_gate, py_gate, py_sim_gate) in zip(
            [CCNot(), CSwap()],
            [braket_ixs.CCNot(controls=pylist([0, 1]), target=2),
             braket_ixs.CSwap(control=0, targets=pylist([1, 2]))],
            [braket_sim_gates.CCNot(targets=pylist([0, 1, 2])),
             braket_sim_gates.CSwap(targets=pylist([0, 1, 2]))],
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1, 2]) == py_gate)
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
        @testset for (jl_noise, py_noise, py_sim_noise) in zip(
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
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, 0)
            @test pyconvert(Braket.Instruction, py_sim_noise) == Braket.Instruction(jl_noise, 0)
        end
        @testset for (jl_noise, py_noise, py_sim_noise) in zip(
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
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, [0, 1])
            @test pyconvert(Braket.Instruction, py_sim_noise) == Braket.Instruction(jl_noise, [0, 1])
        end
    end
    @testset "Result types" begin
        jl_tp_1 = ["h", "x"]
        py_tp_1 = pylist([pystr("h"), pystr("x")])
        jl_herm = Braket.complex_matrix_to_ir(jl_mat)
        py_herm = py_mat
        jl_tp_2 = Union{String, Vector{Vector{Vector{Float64}}}}["h", jl_herm]
        py_tp_2 = pylist([pystr("h"), py_herm])
        @testset for (jl_rt, py_rt) in zip(
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
        )
            @test pyconvert(Bool, Py(jl_rt) == py_rt)
            @test pyconvert(Braket.IR.AbstractProgramResult, py_rt) == jl_rt
        end
    end
    @testset "Result types from default sim" begin
        jl_tp_1 = Braket.Observables.TensorProduct(["h", "x"])
        py_tp_1 = braket_sim_obs.TensorProduct(pylist([braket_sim_obs.Hadamard(targets=pylist([0])), braket_sim_obs.PauliX(targets=pylist([1]))]))
        jl_herm = Braket.Observables.HermitianObservable(jl_mat)
        py_herm = braket_sim_obs.Hermitian(targets=pylist([0]), matrix=py_sim_mat)
        jl_tp_2 = Braket.Observables.TensorProduct([Braket.Observables.H(), jl_herm])
        py_tp_2 = braket_sim_obs.TensorProduct(pylist([braket_sim_obs.Hadamard(targets=pylist([0])), braket_sim_obs.Hermitian(targets=pylist([1]), matrix=py_sim_mat)]))
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
            @test pyconvert(Braket.Result, py_rt) == jl_rt
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
            @test pyconvert(Braket.Result, py_rt) == jl_rt
        end
    end
end
