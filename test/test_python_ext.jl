using Test, PythonCall, BraketSimulator
using PythonCall: pyconvert

@testset "Python integration" begin
    braket_rts       = pyimport("braket.ir.jaqcd.results")
    np               = pyimport("numpy")
    jl_mat           = ComplexF64[0.0 1.0; 1.0 0.0]
    py_mat           = pylist([pylist([pylist([0.0; 0.0]); pylist([1.0; 0.0])]); pylist([pylist([1.0; 0.0]); pylist([0.0; 0.0])])])
    py_sim_mat       = np.array(pylist([pylist([0.0, 1.0]); pylist([1.0, 0.0])]))
    @testset "Result types" begin
        jl_tp_1 = Union{String, Vector{Vector{Vector{Float64}}}}["h", "x"]
        py_tp_1 = pylist([pystr("h"), pystr("x")])
        jl_herm = BraketSimulator.complex_matrix_to_ir(jl_mat)
        py_herm = py_mat
        jl_tp_2 = Union{String, Vector{Vector{Vector{Float64}}}}["h", jl_herm]
        py_tp_2 = pylist([pystr("h"), py_herm])
        @testset for (jl_rt, py_rt, jl_braket_rt) in zip(
                [
                    BraketSimulator.IR.StateVector("statevector"),
                    BraketSimulator.IR.Amplitude(["00", "11"], "amplitude"),
                    BraketSimulator.IR.Probability(nothing, "probability"),
                    BraketSimulator.IR.Probability([0], "probability"),
                    BraketSimulator.IR.DensityMatrix(nothing, "densitymatrix"),
                    BraketSimulator.IR.DensityMatrix([0], "densitymatrix"),
                    BraketSimulator.IR.Sample("x", nothing, "sample"),
                    BraketSimulator.IR.Sample("x", [0], "sample"),
                    BraketSimulator.IR.Sample(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], nothing, "sample"),
                    BraketSimulator.IR.Sample(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], [0], "sample"),
                    BraketSimulator.IR.Sample(jl_tp_1, [0, 1], "sample"),
                    BraketSimulator.IR.Sample(jl_tp_2, [0, 1], "sample"),
                    BraketSimulator.IR.Expectation("x", nothing, "expectation"),
                    BraketSimulator.IR.Expectation("x", [0], "expectation"),
                    BraketSimulator.IR.Expectation(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], nothing, "expectation"),
                    BraketSimulator.IR.Expectation(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], [0], "expectation"),
                    BraketSimulator.IR.Expectation(jl_tp_1, [0, 1], "expectation"),
                    BraketSimulator.IR.Expectation(jl_tp_2, [0, 1], "expectation"),
                    BraketSimulator.IR.Variance("x", nothing, "variance"),
                    BraketSimulator.IR.Variance("x", [0], "variance"),
                    BraketSimulator.IR.Variance(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], nothing, "variance"),
                    BraketSimulator.IR.Variance(Union{Vector{Vector{Vector{Float64}}}, String}[jl_herm], [0], "variance"),
                    BraketSimulator.IR.Variance(jl_tp_1, [0, 1], "variance"),
                    BraketSimulator.IR.Variance(jl_tp_2, [0, 1], "variance"),
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
                    BraketSimulator.StateVector(),
                    BraketSimulator.Amplitude(["00", "11"]),
                    BraketSimulator.Probability(),
                    BraketSimulator.Probability(0),
                    BraketSimulator.DensityMatrix(),
                    BraketSimulator.DensityMatrix(0),
                    BraketSimulator.Sample(BraketSimulator.Observables.X(), Int[]),
                    BraketSimulator.Sample(BraketSimulator.Observables.X(), [0]),
                    BraketSimulator.Sample(BraketSimulator.Observables.HermitianObservable(jl_mat), Int[]),
                    BraketSimulator.Sample(BraketSimulator.Observables.HermitianObservable(jl_mat), [0]),
                    BraketSimulator.Sample(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.X()]), [0, 1]),
                    BraketSimulator.Sample(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.HermitianObservable(jl_mat)]), [0, 1]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.X(), Int[]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.X(), [0]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.HermitianObservable(jl_mat), Int[]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.HermitianObservable(jl_mat), [0]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.X()]), [0, 1]),
                    BraketSimulator.Expectation(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.HermitianObservable(jl_mat)]), [0, 1]),
                    BraketSimulator.Variance(BraketSimulator.Observables.X(), Int[]),
                    BraketSimulator.Variance(BraketSimulator.Observables.X(), [0]),
                    BraketSimulator.Variance(BraketSimulator.Observables.HermitianObservable(jl_mat), Int[]),
                    BraketSimulator.Variance(BraketSimulator.Observables.HermitianObservable(jl_mat), [0]),
                    BraketSimulator.Variance(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.X()]), [0, 1]),
                    BraketSimulator.Variance(BraketSimulator.Observables.TensorProduct([BraketSimulator.Observables.H(), BraketSimulator.Observables.HermitianObservable(jl_mat)]), [0, 1]),
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
        oq3_program = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), sv_adder, Dict("a_in"=>3, "b_in"=>7))
        @test pyconvert(BraketSimulator.OpenQasmProgram, Py(oq3_program)) == oq3_program
        @testset "Full Python circuit execution" begin
            @testset "OpenQASM3" begin
                n_qubits = 5
                # test a "batch"
                sv_simulator = StateVectorSimulator(n_qubits, 0)
                py_inputs    = PyList{Any}([pydict(Dict("a_in"=>2, "b_in"=>5)), pydict(Dict("a_in"=>3, "b_in"=>2))])
                oq3_results  = BraketSimulator.simulate(sv_simulator, PyList{Any}([oq3_program, oq3_program]); inputs=py_inputs, shots=0)
                for oq3_result in oq3_results
                    @test pyconvert(Vector{Float64}, oq3_result.resultTypes[0].value) ≠ pyconvert(Vector{Float64}, oq3_result.resultTypes[1].value)
                end
                # test a "batch" of length 1
                sv_simulator = StateVectorSimulator(n_qubits, 0)
                py_inputs    = PyList{Any}([pydict(Dict("a_in"=>2, "b_in"=>5))])
                oq3_results  = BraketSimulator.simulate(sv_simulator, PyList{Any}([oq3_program]); inputs=py_inputs, shots=0)
                for oq3_result in oq3_results
                    @test pyconvert(Vector{Float64}, oq3_result.resultTypes[0].value) ≠ pyconvert(Vector{Float64}, oq3_result.resultTypes[1].value)
                end
                simple_bell_qasm = """
                h \$0;
                cnot \$0, \$1;
                #pragma braket result state_vector
                """
                sv_simulator = StateVectorSimulator(2, 0)
                oq3_program = BraketSimulator.OpenQasmProgram(BraketSimulator.braketSchemaHeader("braket.ir.openqasm.program", "1"), simple_bell_qasm, nothing)
                oq3_results  = BraketSimulator.simulate(sv_simulator, PyList{Any}([oq3_program]); shots=0)
                for oq3_result in oq3_results
                    @test pyconvert(Vector{ComplexF64}, oq3_result.resultTypes[0].value) ≈ 1/√2 * [1; 0; 0; 1]
                end
            end
        end
    end
end
