using Test, PythonCall, BraketSimulator, Braket

@testset "Python integration" begin
    braket_ixs = pyimport("braket.ir.jaqcd.instructions")
    braket_rts = pyimport("braket.ir.jaqcd.results")
    jl_mat     = ComplexF64[0.0 1.0; 1.0 0.0]
    py_mat     = pylist([pylist([pylist([0.0; 0.0]); pylist([1.0; 0.0])]); pylist([pylist([1.0; 0.0]); pylist([0.0; 0.0])])])
    @testset "Gates" begin
        @testset for (jl_gate, py_gate) in zip(
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
        )
            @test pyconvert(Bool, Py(jl_gate, 0) == py_gate)
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
        end
        angle    = π / 3.5
        angle2   = π / 5.2
        angle3   = π / 0.6
        @testset for (jl_gate, py_gate) in zip(
            [Rx(angle), Ry(angle), Rz(angle), PhaseShift(angle)],
            [
                braket_ixs.Rx(target=0, angle=angle),
                braket_ixs.Ry(target=0, angle=angle),
                braket_ixs.Rz(target=0, angle=angle),
                braket_ixs.PhaseShift(target=0, angle=angle),
            ],
        )
            @test pyconvert(Bool, Py(jl_gate, 0) == py_gate)
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, 0)
        end
        @testset for (jl_gate, py_gate) in zip(
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
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1]) == py_gate)
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset for (jl_gate, py_gate) in zip(
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
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1]) == py_gate)
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1])
        end
        @testset for (jl_gate, py_gate) in zip(
            [CCNot(), CSwap()],
            [braket_ixs.CCNot(controls=pylist([0, 1]), target=2),
             braket_ixs.CSwap(control=0, targets=pylist([1, 2]))],
        )
            @test pyconvert(Bool, Py(jl_gate, [0, 1, 2]) == py_gate)
            @test pyconvert(Braket.Instruction, py_gate) == Braket.Instruction(jl_gate, [0, 1, 2])
        end
    end
    prob = 0.15
    gamma = 0.25
    probx = 0.1
    proby = 0.2
    probz = 0.3
    @testset "Noises" begin
        @testset for (jl_noise, py_noise) in zip(
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
            )
            @test pyconvert(Bool, Py(jl_noise, 0) == py_noise)
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, 0)
        end
        @testset for (jl_noise, py_noise) in zip(
            [TwoQubitDepolarizing(prob), TwoQubitDephasing(prob)],
            [
             braket_ixs.TwoQubitDepolarizing(targets=pylist([0, 1]), probability=prob),
             braket_ixs.TwoQubitDephasing(targets=pylist([0, 1]), probability=prob),
            ],
        )
            @test pyconvert(Bool, Py(jl_noise, [0, 1]) == py_noise)
            @test pyconvert(Braket.Instruction, py_noise) == Braket.Instruction(jl_noise, [0, 1])
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
end
