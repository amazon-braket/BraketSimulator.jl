using BraketSimulator.Quasar, BraketSimulator.Dates, Statistics, LinearAlgebra, BraketSimulator, BraketSimulator.Observables

get_tol(shots::Int) = return (
    shots > 0 ? Dict("atol" => 0.1, "rtol" => 0.15) : Dict("atol" => 0.01, "rtol" => 0)
)
@testset "Quasar" begin
    @testset "Adder" begin
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
        circ = BraketSimulator.to_circuit(sv_adder_qasm, Dict("a_in"=>3, "b_in"=>7))
        correct_instructions = [ 
            BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(6))
            BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(3))
            BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(7))
            BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(4))
            BraketSimulator.Instruction(BraketSimulator.X(), BraketSimulator.QubitSet(8))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 8))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 0))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(0, 8, 4))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 7))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 4))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(4, 7, 3))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 6))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 3))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(3, 6, 2))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 5))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 2))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(2, 5, 1))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 9))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(2, 5, 1))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(1, 2))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 5))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(3, 6, 2))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(2, 3))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 6))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(4, 7, 3))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(3, 4))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 7))
            BraketSimulator.Instruction(BraketSimulator.CCNot(), BraketSimulator.QubitSet(0, 8, 4))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(4, 0))
            BraketSimulator.Instruction(BraketSimulator.CNot(), BraketSimulator.QubitSet(0, 8))
        ]
        for (ix, c_ix) in zip(circ.instructions, correct_instructions)
            @test ix == c_ix
        end
        @test circ.result_types == BraketSimulator.Result[BraketSimulator.Probability(BraketSimulator.QubitSet(9, 5, 6, 7, 8)),
                                                          BraketSimulator.Probability(BraketSimulator.QubitSet(9)),
                                                          BraketSimulator.Probability(BraketSimulator.QubitSet(5, 6, 7, 8))]
    end
    @testset "Measure with density matrix" begin
        qasm = """
        qubit[3] qs;
        qubit q;

        x qs[{0, 2}];
        h q;
        cphaseshift(1) qs, q;
        phaseshift(-2) q;
        measure qs;
        measure q;
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        simulation = BraketSimulator.DensityMatrixSimulator(4, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test BraketSimulator.probabilities(simulation) ≈ (1/2)*[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
    end
    @testset "GPhase" begin
        qasm = """
        qubit[2] qs;

        int[8] two = 2;

        gate x a { U(π, 0, π) a; }
        gate cx c, a { ctrl @ x c, a; }
        gate phase c, a {
            gphase(π/2);
            ctrl(two) @ gphase(π) c, a;
        }
        gate h a { U(π/2, 0, π) a; }

        h qs[0];
        
        cx qs[0], qs[1];
        phase qs[0], qs[1];
        
        gphase(π);
        inv @ gphase(π / 2);
        negctrl @ ctrl @ gphase(2 * π) qs[0], qs[1];
        #pragma braket result amplitude '00', '01', '10', '11'
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        sv = 1/√2 * [-1; 0; 0; 1]
        @test simulation.state_vector ≈ sv 
        @test circuit.result_types == [BraketSimulator.Amplitude(["00", "01", "10", "11"])]
    end
    @testset "Gate def with argument manipulation" begin
        qasm = """
        qubit[2] __qubits__;
        gate u3(θ, ϕ, λ) q {
            gphase(-(ϕ+λ)/2);
            U(θ, ϕ, λ) q;
        }
        u3(0.1, 0.2, 0.3) __qubits__[0];
        """
        circ = BraketSimulator.to_circuit(qasm)
        canonical_ixs = [BraketSimulator.Instruction(BraketSimulator.GPhase{1}(-(0.2 + 0.3)/2), 0), BraketSimulator.Instruction(BraketSimulator.U(0.1, 0.2, 0.3), 0)]
        @test circ.instructions == canonical_ixs 
    end
    @testset "Physical qubits" begin
        qasm = """
        h \$0;
        cnot \$0, \$1;
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        circ = BraketSimulator.to_circuit(visitor)
        @test BraketSimulator.qubit_count(circ) == 2
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.H(), 0), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1])]
    end
    @testset "For loop and subroutines" begin
        qasm_str = """
        OPENQASM 3.0;
        def bell(qubit q0, qubit q1) {
            h q0;
            cnot q0, q1;
        }
        def n_bells(int[32] n, qubit q0, qubit q1) {
            for int i in [0:n - 1] {
                h q0;
                cnot q0, q1;
            }
        }
        qubit[4] __qubits__;
        bell(__qubits__[0], __qubits__[1]);
        n_bells(5, __qubits__[2], __qubits__[3]);
        bit[4] __bit_0__ = "0000";
        __bit_0__[0] = measure __qubits__[0];
        __bit_0__[1] = measure __qubits__[1];
        __bit_0__[2] = measure __qubits__[2];
        __bit_0__[3] = measure __qubits__[3];
        """
        parsed_circ = BraketSimulator.to_circuit(qasm_str, Dict("theta"=>0.2))
        deleteat!(parsed_circ.instructions, length(parsed_circ.instructions)-3:length(parsed_circ.instructions))
        @test parsed_circ.instructions == [BraketSimulator.Instruction(BraketSimulator.H(), 0),
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]),
                                           BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
                                           BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
                                           BraketSimulator.Instruction(BraketSimulator.H(), 2), 
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
                                           BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3]),
                                           BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                           BraketSimulator.Instruction(BraketSimulator.CNot(), [2, 3])
                                          ]
    end
    @testset "Builtin functions" begin
        @testset "Symbolic" begin
            qasm = """
                input float x;
                input float y;
                rx(x) \$0;
                rx(arccos(x)) \$0;
                rx(arcsin(x)) \$0;
                rx(arctan(x)) \$0; 
                rx(ceiling(x)) \$0;
                rx(cos(x)) \$0;
                rx(exp(x)) \$0;
                rx(floor(x)) \$0;
                rx(log(x)) \$0;
                rx(mod(x, y)) \$0;
                rx(sin(x)) \$0;
                rx(sqrt(x)) \$0;
                rx(tan(x)) \$0;
                """
            x = 1.0
            y = 2.0
            inputs      = Dict("x"=>x, "y"=>y)
            parsed_circ = BraketSimulator.to_circuit(qasm, inputs) 
            ixs = [BraketSimulator.Instruction(BraketSimulator.Rx(x), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(acos(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(asin(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(atan(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(ceil(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(cos(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(exp(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(floor(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(log(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(mod(x, y)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(sin(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(sqrt(x)), 0),  
                   BraketSimulator.Instruction(BraketSimulator.Rx(tan(x)), 0)]
            c = BraketSimulator.Circuit()
            for ix in ixs
                BraketSimulator.add_instruction!(c, ix)
            end
            @test parsed_circ.instructions == c.instructions
        end
    end
    @testset "Bad pragma" begin
        qasm = """
        qubit[4] q;
        #pragma braket fake_pragma
        """
        @test_throws Quasar.QasmParseError parse_qasm(qasm)
    end
    @testset "Delay $duration" for (duration, ix) in (("200us", Microsecond(200)),
                                                      ("100μs", Microsecond(100)),
                                                      ("50ns", Nanosecond(50)),
                                                      ("20dt", Nanosecond(20)),
                                                      ("10ms", Millisecond(10)),
                                                      ("1s", Second(1)),
                                                     )
        

        qasm = """
        qubit[4] q;
        x q[0];
        delay[$duration] q[0], q[1];
        """
        @test_warn "delay expression encountered -- currently `delay` is a no-op" parse_qasm(qasm)
        circuit = BraketSimulator.to_circuit(qasm)
        @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.Delay(ix), 0),
                                       BraketSimulator.Instruction(BraketSimulator.Delay(ix), 1),
                                      ]
        simulation = BraketSimulator.StateVectorSimulator(4, 0)
        ref_circ   = BraketSimulator.Circuit()
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.X(), 0))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 1))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 2))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 3))
        ref_sim    = BraketSimulator.StateVectorSimulator(4, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        BraketSimulator.evolve!(ref_sim, ref_circ.instructions)
        @test simulation.state_vector ≈ ref_sim.state_vector
    end
    @testset "Barrier" begin
        qasm = """
        qubit[4] q;
        x q[0];
        barrier q[0];
        #pragma braket result state_vector
        """
        @test_warn "barrier expression encountered -- currently `barrier` is a no-op" parse_qasm(qasm)
        circuit = BraketSimulator.to_circuit(qasm)
        @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Barrier(), 0)]
        simulation = BraketSimulator.StateVectorSimulator(4, 0)
        ref_circ   = BraketSimulator.Circuit()
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.X(), 0))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 1))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 2))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 3))
        ref_sim    = BraketSimulator.StateVectorSimulator(4, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        BraketSimulator.evolve!(ref_sim, ref_circ.instructions)
        @test simulation.state_vector ≈ ref_sim.state_vector
    end
    @testset "Reset" begin
        qasm = """
        qubit[4] q;
        x q[0];
        reset q[0];
        """
        @test_warn "reset expression encountered -- currently `reset` is a no-op" parse_qasm(qasm)
        circuit = BraketSimulator.to_circuit(qasm)
        @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Reset(), 0)]
        simulation = BraketSimulator.StateVectorSimulator(4, 0)
        ref_circ   = BraketSimulator.Circuit()
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.X(), 0))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 1))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 2))
        push!(ref_circ.instructions, BraketSimulator.Instruction(BraketSimulator.I(), 3))
        ref_sim    = BraketSimulator.StateVectorSimulator(4, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        BraketSimulator.evolve!(ref_sim, ref_circ.instructions)
        @test simulation.state_vector ≈ ref_sim.state_vector
        @testset "Reset after Measure" begin
            good_qasm = """
            qubit[4] q;
            measure q[0];
            reset q[0];
            x q[0];
            """
            circuit = BraketSimulator.to_circuit(good_qasm) # no error
            @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.Measure(), 0), BraketSimulator.Instruction(BraketSimulator.Reset(), 0), BraketSimulator.Instruction(BraketSimulator.X(), 0)]
            bad_qasm = """
            qubit[4] q;
            measure q[0];
            reset q[0];
            measure q[0];
            x q[0];
            """
            @test_throws ErrorException("cannot apply instruction to measured qubits.") BraketSimulator.to_circuit(bad_qasm)
        end
    end
    @testset "Adjoint Gradient pragma" begin
        qasm = """
        input float theta;
        qubit[4] q;
        rx(theta) q[0];
        #pragma braket result adjoint_gradient expectation(-6 * y(q[0]) @ i(q[1]) + 0.75 * y(q[2]) @ z(q[3])) theta
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor(Dict("theta"=>0.1))
        @test_throws Quasar.QasmVisitorError("Result type adjoint_gradient is not supported.", "TypeError") visitor(parsed)
        #= circuit = BraketSimulator.Circuit(qasm, Dict("theta"=>0.1))
        θ       = FreeParameter("theta")
        obs     = -2 * BraketSimulator.Observables.Y() * (3 * BraketSimulator.Observables.I()) + 0.75 * BraketSimulator.Observables.Y() * BraketSimulator.Observables.Z()
        @test circuit.result_types[1].observable == obs
        @test circuit.result_types[1].targets == [QubitSet([0, 1]), QubitSet([2, 3])]
        @test circuit.result_types[1].parameters == ["theta"]=#
    end
    @testset "Switch/case" begin
        qasm = """
        input int[8] x;
        switch (x + 1) {
            case 0b00 {}
            default { z \$0; }
        }
        """
        circ = BraketSimulator.to_circuit(qasm, Dict("x"=> -1))
        @test isempty(circ.instructions)
        circ = BraketSimulator.to_circuit(qasm, Dict("x"=> 0))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        qasm = """
        input int[8] x;
        switch (x) { case 0 {} case 1, 2 { z \$0; }  }
        """
        circ = BraketSimulator.to_circuit(qasm, Dict("x"=> 0))
        @test isempty(circ.instructions)
        circ = BraketSimulator.to_circuit(qasm, Dict("x"=> 1))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        circ = BraketSimulator.to_circuit(qasm, Dict("x"=> 2))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        qasm = """
        input int[8] x;
        switch (x) { case 0 {} default { z \$0; } default { x \$0; } }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.to_circuit(qasm, Dict("x"=>0))
        qasm = """
        input int[8] x;
        switch (x) { default { z \$0; } case 0 {} }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.to_circuit(qasm, Dict("x"=>0))
        qasm = """
        input int[8] x;
        switch (x) { case 0 { z \$0; } true {} }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.to_circuit(qasm, Dict("x"=>0))
    end
    @testset "Global gate control" begin
        qasm = """
        qubit q1;
        qubit q2;

        h q1;
        h q2;
        ctrl @ s q1, q2;
        """
        circuit    = BraketSimulator.to_circuit(qasm) 
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ [0.5, 0.5, 0.5, 0.5im]
    end
    @testset "Simple Pow" for (orig_gate, sqrt_gate, sqrt_ix) in (("z", "s", BraketSimulator.S()),
                                                                  ("x", "v", BraketSimulator.V()),
                                                                 )
        custom_qasm = """
        qubit q1;
        qubit q2;
        h q1;
        h q2;
        
        pow(1/2) @ $orig_gate q1;
        """
        standard_qasm = """
        qubit q1;
        qubit q2;
        h q1;
        h q2;
        
        $sqrt_gate q1;
        """
        canonical_ixs = [BraketSimulator.Instruction(BraketSimulator.H(), 0),
                         BraketSimulator.Instruction(BraketSimulator.H(), 1),
                         BraketSimulator.Instruction(sqrt_ix, 0)]
        canonical_simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(canonical_simulation, canonical_ixs)
        @testset "$title" for (title, qasm) in (("standard", standard_qasm), ("custom", custom_qasm))
            circuit    = BraketSimulator.to_circuit(qasm)
            simulation = BraketSimulator.StateVectorSimulator(2, 0)
            BraketSimulator.evolve!(simulation, circuit.instructions)
            @test simulation.state_vector ≈ canonical_simulation.state_vector
        end
    end
    @testset "Complex Pow" begin
        custom_qasm = """
        int[8] two = 2;
        gate x a { U(π, 0, π) a; }
        gate cx c, a {
            pow(1) @ ctrl @ x c, a;
        }
        gate cxx_1 c, a {
            pow(two) @ cx c, a;
        }
        gate cxx_2 c, a {
            pow(1/2) @ pow(4) @ cx c, a;
        }
        gate cxxx c, a {
            pow(1) @ pow(two) @ cx c, a;
        }

        qubit q1;
        qubit q2;
        qubit q3;
        qubit q4;
        qubit q5;

        pow(1/2) @ x q1;       // half flip
        pow(1/2) @ x q1;       // half flip
        cx q1, q2;   // flip
        cxx_1 q1, q3;    // don't flip
        cxx_2 q1, q4;    // don't flip
        cnot q1, q5;    // flip
        x q3;       // flip
        x q4;       // flip

        s q1;   // sqrt z
        s q1;   // again
        inv @ z q1; // inv z
        """
        standard_qasm = """
        int[8] two = 2;
        gate cxx_1 c, a {
            pow(two) @ cnot c, a;
        }
        gate cxx_2 c, a {
            pow(1/2) @ pow(4) @ cnot c, a;
        }
        gate cxxx c, a {
            pow(1) @ pow(two) @ cnot c, a;
        }

        qubit q1;
        qubit q2;
        qubit q3;
        qubit q4;
        qubit q5;

        pow(1/2) @ x q1;       // half flip
        pow(1/2) @ x q1;       // half flip
        cnot q1, q2;   // flip
        cxx_1 q1, q3;    // don't flip
        cxx_2 q1, q4;    // don't flip
        cnot q1, q5;    // flip
        x q3;       // flip
        x q4;       // flip

        s q1;   // sqrt z
        s q1;   // again
        inv @ z q1; // inv z
        """
        @testset "$title" for (title, qasm) in (("standard", standard_qasm),
                                                ("custom", custom_qasm)
                                               )
            circuit    = BraketSimulator.to_circuit(qasm) 
            simulation = BraketSimulator.StateVectorSimulator(5, 0)
            BraketSimulator.evolve!(simulation, circuit.instructions)
            sv = zeros(32)
            sv[end] = 1.0
            @test simulation.state_vector ≈ sv
        end
    end
    @testset "Gate control" begin
        qasm = """
        int[8] two = 2;
        gate x a { U(π, 0, π) a; }
        gate cx c, a {
            ctrl @ x c, a;
        }
        gate ccx_1 c1, c2, a {
            ctrl @ ctrl @ x c1, c2, a;
        }
        gate ccx_2 c1, c2, a {
            ctrl(two) @ x c1, c2, a;
        }
        gate ccx_3 c1, c2, a {
            ctrl @ cx c1, c2, a;
        }

        qubit q1;
        qubit q2;
        qubit q3;
        qubit q4;
        qubit q5;

        // doesn't flip q2
        cx q1, q2;
        // flip q1
        x q1;
        // flip q2
        cx q1, q2;
        // doesn't flip q3, q4, q5
        ccx_1 q1, q4, q3;
        ccx_2 q1, q3, q4;
        ccx_3 q1, q3, q5;
        // flip q3, q4, q5;
        ccx_1 q1, q2, q3;
        ccx_2 q1, q2, q4;
        ccx_2 q1, q2, q5;
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        expected_ixs = [
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,)), [0, 1]),
                        BraketSimulator.Instruction(BraketSimulator.U(π, 0, π), 0),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,)), [0, 1]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 3, 2]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 2, 3]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 2, 4]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 1, 2]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 1, 3]),
                        BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.U(π, 0, π), (1,1)), [0, 1, 4]),
                       ]
        @test length(circuit.instructions) == length(expected_ixs)
        @testset for ix in 1:length(expected_ixs)
            @test circuit.instructions[ix] == expected_ixs[ix]
        end
        simulation = BraketSimulator.StateVectorSimulator(5, 0)
        sv = zeros(ComplexF64, 32)
        sv[end] = 1.0
        canonical_ixs        = [
                                BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]),
                                BraketSimulator.Instruction(BraketSimulator.X(), 0),
                                BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 3, 2]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 2, 3]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 2, 4]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 1, 2]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 1, 3]),
                                BraketSimulator.Instruction(BraketSimulator.CCNot(), [0, 1, 4])
                               ]
        canonical_simulation = BraketSimulator.StateVectorSimulator(5, 0)
        ii = 1
        @testset for (ix, c_ix) in zip(circuit.instructions, canonical_ixs)
            BraketSimulator.evolve!(simulation, [ix])
            BraketSimulator.evolve!(canonical_simulation, [c_ix])
            for jj in 1:32
                @test simulation.state_vector[jj] ≈ canonical_simulation.state_vector[jj] atol=1e-10
            end
            @test simulation.state_vector ≈ canonical_simulation.state_vector atol=1e-10
            ii += 1
        end
        @test canonical_simulation.state_vector ≈ sv rtol=1e-10
        @test simulation.state_vector ≈ sv rtol=1e-10
        qasm = """
        gate cccx c1, c2, c3, a {
            ctrl(3) @ x c1, c2, c3, a;
        }
        gate ncccx c1, c2, c3, a {
            negctrl(3) @ x c1, c2, c3, a;
        }

        qubit q1;
        qubit q2;
        qubit q3;
        qubit q4;
        qubit q5;
        
        h q1;
        h q2;
        h q3;
        h q4;
        h q5;
        cccx q1, q2, q5, q3;
        ncccx q4, q2, q5, q3;
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.H(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.H(), 1),
                                       BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                       BraketSimulator.Instruction(BraketSimulator.H(), 3),
                                       BraketSimulator.Instruction(BraketSimulator.H(), 4),
                                       BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.X(), (1, 1, 1)), [0, 1, 4, 2]),
                                       BraketSimulator.Instruction(BraketSimulator.Control(BraketSimulator.X(), (0, 0, 0)), [3, 1, 4, 2]),
                                      ]
    end
    @testset "Gate inverses" begin
        qasm = """
        gate rand_u_1 a { U(1, 2, 3) a; }
        gate rand_u_2 a { U(2, 3, 4) a; }
        gate rand_u_3 a { inv @ U(3, 4, 5) a; }

        gate both a {
            rand_u_1 a;
            rand_u_2 a;
        }
        gate both_inv a {
            inv @ both a;
        }
        gate all_3 a {
            rand_u_1 a;
            rand_u_2 a;
            rand_u_3 a;
        }
        gate all_3_inv a {
            inv @ inv @ inv @ all_3 a;
        }

        gate apply_phase a {
            gphase(1);
        }

        gate apply_phase_inv a {
            inv @ gphase(1);
        }

        qubit q;

        both q;
        both_inv q;

        all_3 q;
        all_3_inv q;

        apply_phase q;
        apply_phase_inv q;

        U(1, 2, 3) q;
        inv @ U(1, 2, 3) q;

        s q;
        inv @ s q;

        t q;
        inv @ t q;
        """
        circuit   = BraketSimulator.to_circuit(qasm) 
        collapsed = prod(BraketSimulator.matrix_rep(ix.operator) for ix in circuit.instructions)
        @test collapsed ≈ diagm(ones(ComplexF64, 2^BraketSimulator.qubit_count(circuit)))
    end
    @testset "Noise pragmas $noise, $qubits" for (noise, arg, qubits, ix) in (("bit_flip", ".5", "qs[1]", BraketSimulator.Instruction(BraketSimulator.BitFlip(0.5), 1)),
                                                                     ("phase_flip", ".5", "qs[0]", BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.5), 0)),
                                                                     ("pauli_channel", ".1, .2, .3", "qs[0]", BraketSimulator.Instruction(BraketSimulator.PauliChannel(0.1, 0.2, 0.3), 0)),
                                                                     ("depolarizing", ".5", "qs[0]", BraketSimulator.Instruction(BraketSimulator.Depolarizing(0.5), 0)),
                                                                     ("two_qubit_depolarizing", ".9", "qs", BraketSimulator.Instruction(BraketSimulator.TwoQubitDepolarizing(0.9), [0, 1])),
                                                                     ("two_qubit_depolarizing", ".7", "qs[1], qs[0]", BraketSimulator.Instruction(BraketSimulator.TwoQubitDepolarizing(0.7), [1, 0])),
                                                                     ("two_qubit_dephasing", ".6", "qs", BraketSimulator.Instruction(BraketSimulator.TwoQubitDephasing(0.6), [0, 1])),
                                                                     ("amplitude_damping", ".2", "qs[0]", BraketSimulator.Instruction(BraketSimulator.AmplitudeDamping(0.2), 0)),
                                                                     ("generalized_amplitude_damping", ".2, .3", "qs[1]", BraketSimulator.Instruction(BraketSimulator.GeneralizedAmplitudeDamping(0.2, 0.3), 1)),
                                                                     ("phase_damping", ".4", "qs[0]", BraketSimulator.Instruction(BraketSimulator.PhaseDamping(0.4), 0)),
                                                                     ("kraus", "[[0.9486833im, 0], [0, 0.9486833im]], [[0, 0.31622777], [0.31622777, 0]]", "qs[0]", BraketSimulator.Instruction(BraketSimulator.Kraus([[0.9486833im 0; 0 0.9486833im], [0 0.31622777; 0.31622777 0]]), 0)),
                                                                     ("kraus", "[[0.9486832980505138, 0, 0, 0], [0, 0.9486832980505138, 0, 0], [0, 0, 0.9486832980505138, 0], [0, 0, 0, 0.9486832980505138]], [[0, 0.31622776601683794, 0, 0], [0.31622776601683794, 0, 0, 0], [0, 0, 0, 0.31622776601683794], [0, 0, 0.31622776601683794, 0]]", "qs[{1, 0}]", BraketSimulator.Instruction(BraketSimulator.Kraus([[0.9486832980505138 0 0 0; 0 0.9486832980505138 0 0; 0 0 0.9486832980505138 0; 0 0 0 0.9486832980505138], [0 0.31622776601683794 0 0; 0.31622776601683794 0 0 0; 0 0 0 0.31622776601683794; 0 0 0.31622776601683794 0]]), [1, 0])),
                                                                    )

        
        qasm = """
        qubit[2] qs;

        #pragma braket noise $noise($arg) $qubits
        """
        circuit   = BraketSimulator.to_circuit(qasm) 
        @test only(circuit.instructions) == ix
    end
    @testset "StandardObservable/target mismatch" begin
        @testset "Observable $obs, result type $rt, qubits $qubits" for obs in ("h", "x", "y", "z"),
                                                                        rt in ("expectation", "variance", "sample"),
                                                                        qubits in ("q[1:2]", "q[0:1]", "q[{0, 2}]")
            qasm = """
            qubit[3] q;
            i q;
            
            #pragma braket result $rt $obs($qubits)
            """
            parsed  = parse_qasm(qasm)
            visitor = Quasar.QasmProgramVisitor()
            @test_throws Quasar.QasmVisitorError("Standard observable target must be exactly 1 qubit.", "ValueError") visitor(parsed)
        end
    end
    @testset "Basis rotations" begin
        @testset "StandardObservables" begin
            qasm = """
            qubit[3] q;
            i q;
            
            #pragma braket result expectation z(q[2]) @ x(q[0])
            #pragma braket result variance x(q[0]) @ y(q[1])
            #pragma braket result sample x(q[0])
            """
            circuit = BraketSimulator.to_circuit(qasm)
            BraketSimulator.basis_rotation_instructions!(circuit)
            c_bris  = [circuit.basis_rotation_instructions[1], BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(mapreduce(ix->BraketSimulator.matrix_rep(ix.operator), *, circuit.basis_rotation_instructions[2:end]))), [1])]
            bris   = vcat(BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.diagonalizing_gates(BraketSimulator.Observables.Y(), [1]))
            for (ix, bix) in zip(c_bris, bris)
                @test BraketSimulator.matrix_rep(ix.operator) ≈ transpose(BraketSimulator.matrix_rep(bix.operator))
                @test ix.target == bix.target
            end
        end
        @testset "Identity" begin
            qasm = """
            qubit[3] q;
            i q;
            
            #pragma braket result expectation z(q[2]) @ x(q[0])
            #pragma braket result variance x(q[0]) @ y(q[1])
            #pragma braket result sample i(q[0])
            """
            circuit = BraketSimulator.to_circuit(qasm) 
            BraketSimulator.basis_rotation_instructions!(circuit)
            c_bris  = [circuit.basis_rotation_instructions[1], BraketSimulator.Instruction(BraketSimulator.Unitary(Matrix(mapreduce(ix->BraketSimulator.matrix_rep(ix.operator), *, circuit.basis_rotation_instructions[2:end]))), [1])]
            bris   = vcat(BraketSimulator.Instruction(BraketSimulator.H(), [0]), BraketSimulator.diagonalizing_gates(BraketSimulator.Observables.Y(), [1]))
            for (ix, bix) in zip(c_bris, bris)
                @test BraketSimulator.matrix_rep(ix.operator) ≈ transpose(BraketSimulator.matrix_rep(bix.operator))
                @test ix.target == bix.target
            end
        end
        @testset "Hermitian" begin
            qasm = """
            qubit[3] q;
            i q;
            #pragma braket result expectation x(q[2])
            // # noqa: E501
            #pragma braket result expectation hermitian([[-6+0im, 2+1im, -3+0im, -5+2im], [2-1im, 0im, 2-1im, -5+4im], [-3+0im, 2+1im, 0im, -4+3im], [-5-2im, -5-4im, -4-3im, -6+0im]]) q[0:1]
            // # noqa: E501
            #pragma braket result expectation x(q[2]) @ hermitian([[-6+0im, 2+1im, -3+0im, -5+2im], [2-1im, 0im, 2-1im, -5+4im], [-3+0im, 2+1im, 0im, -4+3im], [-5-2im, -5-4im, -4-3im, -6+0im]]) q[0:1]
            """
            circuit = BraketSimulator.to_circuit(qasm) 
            BraketSimulator.basis_rotation_instructions!(circuit)
            arr = [-6 2+1im -3 -5+2im;
                    2-1im 0 2-1im -5+4im;
                   -3 2+1im 0 -4+3im;
                   -5-2im -5-4im -4-3im -6]
            h = BraketSimulator.Observables.HermitianObservable(arr)
            bris = vcat(BraketSimulator.diagonalizing_gates(h, [0, 1]), BraketSimulator.Instruction(BraketSimulator.H(), [2]))
            for (ix, bix) in zip(circuit.basis_rotation_instructions, bris)
                @test Matrix(BraketSimulator.matrix_rep(ix.operator)) ≈ adjoint(BraketSimulator.fix_endianness(Matrix(BraketSimulator.matrix_rep(bix.operator))))
                @test ix.target == bix.target
            end
        end
    end
    @testset "Unitary pragma" begin
        standard_qasm = """
        qubit[3] q;

        x q[0];
        h q[1];

        t q[0];
        ti q[0];

        h q[1];
        h q[1];

        ccnot q;
        """
        unitary_qasm = """
        qubit[3] q;

        x q[0];
        h q[1];

        // unitary pragma for t gate
        #pragma braket unitary([[1.0, 0], [0, 0.70710678 + 0.70710678im]]) q[0]
        ti q[0];

        // unitary pragma for h gate (with phase shift)
        #pragma braket unitary([[0.70710678im, 0.70710678im], [0 - -0.70710678im, -0.0 - 0.70710678im]]) q[1]
        gphase(-π/2) q[1];
        h q[1];

        // unitary pragma for ccnot gate
        #pragma braket unitary([[1.0, 0, 0, 0, 0, 0, 0, 0], [0, 1.0, 0, 0, 0, 0, 0, 0], [0, 0, 1.0, 0, 0, 0, 0, 0], [0, 0, 0, 1.0, 0, 0, 0, 0], [0, 0, 0, 0, 1.0, 0, 0, 0], [0, 0, 0, 0, 0, 1.0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1.0], [0, 0, 0, 0, 0, 0, 1.0, 0]]) q
        """
        circuit     = BraketSimulator.to_circuit(unitary_qasm)
        ref_circuit = BraketSimulator.to_circuit(standard_qasm)
        simulation  = BraketSimulator.StateVectorSimulator(3, 0)
        ref_sim     = BraketSimulator.StateVectorSimulator(3, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        BraketSimulator.evolve!(ref_sim, ref_circuit.instructions)
        @test simulation.state_vector ≈ ref_sim.state_vector
    end
    @testset "Input" begin
        qasm = """
        input int[8] in_int;
        input bit[8] in_bit;
        int[8] doubled;

        doubled = in_int * 2;
        """
        in_bit = 0b10110010
        parsed  = parse_qasm(qasm)
        @testset for in_int in (0, 1, -2, 5)
            inputs  = Dict("in_int"=>in_int, "in_bit"=>in_bit)
            visitor = QasmProgramVisitor(inputs)
            visitor(parsed)
            @test visitor.classical_defs["doubled"].val == in_int * 2
            @test visitor.classical_defs["in_bit"].val  == in_bit
        end
    end
    @testset "Gate on qubit registers" begin 
        qasm = """
        qubit[3] qs;
        qubit q;

        x qs[{0, 2}];
        h q;
        cphaseshift(1) qs, q;
        phaseshift(-2) q;
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        simulation = BraketSimulator.StateVectorSimulator(4, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ (1/√2)*[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
    end
    @testset "Verbatim" begin
        with_verbatim = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[2] q;
        #pragma braket verbatim
        box{
            cnot q[0], q[1];
            cnot q[0], q[1];
            rx(1.57) q[0];
        }
        b[0] = measure q[0];
        b[1] = measure q[1];
        """
        circuit        = BraketSimulator.to_circuit(with_verbatim)
        sim_w_verbatim = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(sim_w_verbatim, circuit.instructions)

        without_verbatim = """
        OPENQASM 3.0;
        bit[2] b;
        qubit[2] q;
        box{
            cnot q[0], q[1];
            cnot q[0], q[1];
            rx(1.57) q[0];
        }
        b[0] = measure q[0];
        b[1] = measure q[1];
        """
        circuit         = BraketSimulator.to_circuit(without_verbatim)
        sim_wo_verbatim = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(sim_wo_verbatim, circuit.instructions)

        @test sim_w_verbatim.state_vector ≈ sim_wo_verbatim.state_vector
    end
    @testset "Void subroutine" begin
        qasm = """
               def flip(qubit q) {
                 x q;
               }
               qubit[2] qs;
               flip(qs[0]);
               """
        circuit    = BraketSimulator.to_circuit(qasm)
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ [0, 0, 1, 0]
    end
    @testset "Rotation parameter expressions" begin
        @testset  "Operation: $operation" for (operation, state_vector) in
            [
                ["rx(π) q[0];", [0, -im]],
                ["rx(pi) q[0];", [0, -im]],
                ["rx(ℇ) q[0];", [0.21007866, -0.97768449im]],
                ["rx(euler) q[0];", [0.21007866, -0.97768449im]],
                ["rx(τ) q[0];", [-1, 0]],
                ["rx(tau) q[0];", [-1, 0]],
                ["rx(pi + pi) q[0];", [-1, 0]],
                ["rx(pi - pi) q[0];", [1, 0]],
                ["rx(-pi + pi) q[0];", [1, 0]],
                ["rx(pi * 2) q[0];", [-1, 0]],
                ["rx(pi / 2) q[0];", [0.70710678, -0.70710678im]],
                ["rx(-pi / 2) q[0];", [0.70710678, 0.70710678im]],
                ["rx(-pi) q[0];", [0, im]],
                ["rx(pi + 2 * pi) q[0];", [0, im]],
                ["rx(pi + pi / 2) q[0];", [-0.70710678, -0.70710678im]],
                ["rx((pi / 4) + (pi / 2) / 2) q[0];", [0.70710678, -0.70710678im]],
                ["rx(0) q[0];", [1, 0]],
                ["rx(0 + 0) q[0];", [1, 0]],
                ["rx((1.1 + 2.04) / 2) q[0];", [0.70738827, -0.70682518im]],
                ["rx((6 - 2.86) * 0.5) q[0];", [0.70738827, -0.70682518im]],
                ["rx(pi ** 2) q[0];", [0.22058404, 0.97536797im]],
            ]
            qasm = """
            OPENQASM 3.0;
            bit[1] b;
            qubit[1] q;
            $operation
            #pragma braket result state_vector
            """
            circuit    = BraketSimulator.to_circuit(qasm) 
            simulation = BraketSimulator.StateVectorSimulator(1, 0)
            BraketSimulator.evolve!(simulation, circuit.instructions)
            @test simulation.state_vector ≈ state_vector
        end
    end
    @testset "No result types no shots" begin
        qasm = """
        qubit[2] q;

        h q[0];
        ctrl @ x q[0], q[1];
        """
        circuit    = BraketSimulator.to_circuit(qasm)
        program    = BraketSimulator.Program(circuit)
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        @test_throws BraketSimulator.ValidationError BraketSimulator.simulate(simulation, program, 0)
    end
    @testset "Referencing invalid qubit(s)" begin
        source = """
        qubit[2] q;
        z q;
        #pragma braket result expectation z(q[2])
        """
        parsed  = parse_qasm(source)
        visitor = QasmProgramVisitor()
        @test_throws Quasar.QasmVisitorError("Invalid qubit index '2' in 'q'.", "IndexError") visitor(parsed)
        set_source = """
        qubit[4] q;
        z q;
        #pragma braket result expectation z(q[{0,2,5}])
        """
        parsed  = parse_qasm(set_source)
        visitor = QasmProgramVisitor()
        @test_throws Quasar.QasmVisitorError("Invalid qubit index '5' in 'q'.", "IndexError") visitor(parsed)
    end
    @testset "Invalid Hermitian target" begin
        qasm = """
        OPENQASM 3.0;
        qubit[3] q;
        i q;
        #pragma braket result expectation hermitian([[-6+0im, 2+1im, -3+0im, -5+2im], [2-1im, 0im, 2-1im, -5+4im], [-3+0im, 2+1im, 0im, -4+3im], [-5-2im, -5-4im, -4-3im, -6+0im]]) q[0]
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        err_msg = "Invalid observable specified: [[[-6.0, 0.0], [2.0, 1.0], [-3.0, 0.0], [-5.0, 2.0]], [[2.0, -1.0], [0.0, 0.0], [2.0, -1.0], [-5.0, 4.0]], [[-3.0, 0.0], [2.0, 1.0], [0.0, 0.0], [-4.0, 3.0]], [[-5.0, -2.0], [-5.0, -4.0], [-4.0, -3.0], [-6.0, 0.0]]], targets: [0]" 
        @test_throws Quasar.QasmVisitorError(err_msg, "ValueError") visitor(parsed)
    end
    @testset "GRCS 16" begin
        simulator = BraketSimulator.StateVectorSimulator(16, 0)
        circuit   = BraketSimulator.to_circuit(joinpath(@__DIR__, "grcs_16.qasm")) 
        BraketSimulator.evolve!(simulator, circuit.instructions)
        probs     = BraketSimulator.probabilities(simulator)        
        @test first(probs) ≈ 0.0000062 atol=1e-7
    end
end
