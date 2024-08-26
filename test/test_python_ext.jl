using Test, JSON3, PythonCall, BraketSimulator

@testset "Python integration" begin
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
        @testset "Full Python circuit execution" begin
            @testset "OpenQASM3" begin
                n_qubits = 5
                # test a "batch"
                py_inputs    = JSON3.write([Dict{String, Any}("a_in"=>2, "b_in"=>5), Dict{String, Any}("a_in"=>3, "b_in"=>2)])
                oq3_results  = BraketSimulator.simulate("braket_sv_v2", PyList([sv_adder, sv_adder]), py_inputs, 0)
                for oq3_result in oq3_results
                    result = JSON3.read(oq3_result, BraketSimulator.GateModelTaskResult)
                    @test result.additionalMetadata.action.source == sv_adder 
                    @test length(result.resultTypes) == 3
                end
                # test a "batch" of length 1
                py_inputs    = JSON3.write([Dict{String, Any}("a_in"=>2, "b_in"=>5)])
                oq3_results  = BraketSimulator.simulate("braket_sv_v2", PyList([sv_adder]), py_inputs, 0)
                for oq3_result in oq3_results
                    result = JSON3.read(oq3_result, BraketSimulator.GateModelTaskResult)
                    @test result.additionalMetadata.action.source == sv_adder 
                    @test length(result.resultTypes) == 3
                end
                simple_bell_qasm = """
                h \$0;
                cnot \$0, \$1;
                #pragma braket result state_vector
                """
                oq3_results  = BraketSimulator.simulate("braket_sv_v2", PyList([simple_bell_qasm]), "[{}]", 0)
                for oq3_result in oq3_results
                    result = JSON3.read(oq3_result, BraketSimulator.GateModelTaskResult)
                    @test result.additionalMetadata.action.source == simple_bell_qasm
                    @test length(result.resultTypes) == 1
                    @test result.resultTypes[1].type == BraketSimulator.IR.StateVector("statevector")
                end
                simple_bell_qasm = """
                h \$0;
                cy \$0, \$1;
                #pragma braket result amplitude '00', '01', '10', '11' 
                """
                oq3_results  = BraketSimulator.simulate("braket_sv_v2", PyList([simple_bell_qasm]), "[{}]", 0)
                for oq3_result in oq3_results
                    result = JSON3.read(oq3_result, BraketSimulator.GateModelTaskResult)
                    @test result.additionalMetadata.action.source == simple_bell_qasm
                    @test length(result.resultTypes) == 1
                    @test result.resultTypes[1].type.type == "amplitude"
                    @test result.resultTypes[1].type.states == ["00", "01", "10", "11"]
                    @test result.resultTypes[1].value == Dict("00"=>[1/√2, 0.0], "01"=>[0.0, 0.0], "10"=>[0.0, 0.0], "11"=>[0.0, 1/√2])
                end
                simple_bell_qasm = """
                h \$0;
                cy \$0, \$1;
                #pragma braket result expectation z(\$0)
                """
                oq3_results  = BraketSimulator.simulate("braket_sv_v2", PyList([simple_bell_qasm]), "[{}]", 0)
                for oq3_result in oq3_results
                    result = JSON3.read(oq3_result, BraketSimulator.GateModelTaskResult)
                    @test result.additionalMetadata.action.source == simple_bell_qasm
                    @test length(result.resultTypes) == 1
                    @test result.resultTypes[1].type.type    == "expectation"
                    @test result.resultTypes[1].type.targets == [0] 
                    @test result.resultTypes[1].type.observable == ["z"]
                    @test result.resultTypes[1].value == 0.0
                end
            end
        end
    end
end
