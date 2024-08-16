using Test, PythonCall, BraketSimulator, Braket
using PythonCall: pyconvert

@testset "Python integration" begin
    braket_rts       = pyimport("braket.ir.jaqcd.results")
    np               = pyimport("numpy")
    jl_mat           = ComplexF64[0.0 1.0; 1.0 0.0]
    py_mat           = pylist([pylist([pylist([0.0; 0.0]); pylist([1.0; 0.0])]); pylist([pylist([1.0; 0.0]); pylist([0.0; 0.0])])])
    py_sim_mat       = np.array(pylist([pylist([0.0, 1.0]); pylist([1.0, 0.0])]))
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
            @test PythonCall.pyconvert(Bool, Py(jl_rt) == py_rt)
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
        @test pyconvert(Braket.OpenQasmProgram, Py(oq3_program)) == oq3_program
        @testset "Full Python circuit execution" begin
            @testset "OpenQASM3" begin
                # test a "batch"
                sv_simulator = StateVectorSimulator(n_qubits, 0)
                py_inputs    = PyList{Any}([pydict(Dict("a_in"=>2, "b_in"=>5)), pydict(Dict("a_in"=>3, "b_in"=>2))])
                oq3_results  = simulate(sv_simulator, PyList{Any}([oq3_program, oq3_program]); inputs=py_inputs, shots=0)
                for oq3_result in oq3_results
                    @test pyconvert(Vector{Float64}, oq3_result.resultTypes[0].value) ≠ pyconvert(Vector{Float64}, oq3_result.resultTypes[1].value)
                end
                # test a "batch" of length 1
                sv_simulator = StateVectorSimulator(n_qubits, 0)
                py_inputs    = PyList{Any}([pydict(Dict("a_in"=>2, "b_in"=>5))])
                oq3_results  = simulate(sv_simulator, PyList{Any}([oq3_program]); inputs=py_inputs, shots=0)
                for oq3_result in oq3_results
                    @test pyconvert(Vector{Float64}, oq3_result.resultTypes[0].value) ≠ pyconvert(Vector{Float64}, oq3_result.resultTypes[1].value)
                end
            end
        end
    end
end
