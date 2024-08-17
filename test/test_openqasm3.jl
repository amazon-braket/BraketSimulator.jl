using BraketSimulator.Quasar, Statistics, LinearAlgebra, BraketSimulator, BraketSimulator.Observables

get_tol(shots::Int) = return (
    shots > 0 ? Dict("atol" => 0.1, "rtol" => 0.15) : Dict("atol" => 0.01, "rtol" => 0)
)
@testset "Quasar" begin
    @testset "QasmExpression" begin
        @testset "Printing" begin
            expr = Quasar.QasmExpression(:head, Quasar.QasmExpression(:integer_literal, 2))
            @test sprint(show, expr) == "QasmExpression :head\n└─ QasmExpression :integer_literal\n   └─ 2\n"
            bad_expression = Quasar.QasmExpression(:invalid_expression, Quasar.QasmExpression(:error))
            err = Quasar.QasmVisitorError("unable to evaluate expression $bad_expression")
            @test sprint(showerror, err) == "QasmVisitorError: unable to evaluate expression $bad_expression"
        end
        @testset "Expression basics" begin
            expr = Quasar.QasmExpression(:head, Quasar.QasmExpression(:integer_literal, 2))
            @test length(expr) == 1
            push!(expr, Quasar.QasmExpression(:integer_literal, 3))
            @test length(expr.args) == 2
            append!(expr, Quasar.QasmExpression(:integer_literal, 4))
            @test length(expr.args) == 3 
            append!(expr, [Quasar.QasmExpression(:integer_literal, 5), Quasar.QasmExpression(:integer_literal, 6)])
            @test length(expr.args) == 5
        end
        @testset "Equality" begin
            expr_a = Quasar.QasmExpression(:head, Quasar.QasmExpression(:integer_literal, 2))
            expr_b = Quasar.QasmExpression(:head, Quasar.QasmExpression(:float_literal, 2.2))
            @test expr_a != expr_b
            expr_a = Quasar.QasmExpression(:head, Quasar.QasmExpression(:integer_literal, 2))
            @test expr_a == copy(expr_a) 
        end
        @testset "Name undefined" begin
            expr = Quasar.QasmExpression(:undef, Quasar.QasmExpression(:integer_literal, 2))
            @test_throws Quasar.QasmVisitorError("name not defined for expressions of type undef") Quasar.name(expr)
        end
    end
    @testset "Visiting/evaluating an invalid expression" begin
        visitor = Quasar.QasmProgramVisitor()
        bad_expression = Quasar.QasmExpression(:invalid_expression, Quasar.QasmExpression(:error))
        @test_throws Quasar.QasmVisitorError("cannot visit expression $bad_expression.") visitor(bad_expression)
        @test_throws Quasar.QasmVisitorError("unable to evaluate expression $bad_expression") Quasar.evaluate(visitor, bad_expression)
        bad_expression = Quasar.QasmExpression(:invalid_expression, Quasar.QasmExpression(:error))
    end
    @testset "Visitors and parents" begin
        qasm = """
               def flip(qubit q) {
                 x q;
               }
               qubit[2] qs;
               flip(qs[0]);
               """
        parsed = parse_qasm(qasm)
        visitor = Quasar.QasmProgramVisitor()
        visitor(parsed)
        for_visitor = Quasar.QasmForLoopVisitor(visitor)
        @test Quasar.function_defs(for_visitor) == Quasar.function_defs(visitor)
        @test Quasar.qubit_defs(for_visitor) == Quasar.qubit_defs(visitor)
        @test BraketSimulator.qubit_count(for_visitor) == Quasar.qubit_count(visitor)
        push!(for_visitor, BraketSimulator.Amplitude("00"))
        @test visitor.results[1] == BraketSimulator.Amplitude("00")
        push!(for_visitor, [BraketSimulator.StateVector(), BraketSimulator.Probability()])
        @test visitor.results[2] == BraketSimulator.StateVector()
        @test visitor.results[3] == BraketSimulator.Probability()
        push!(for_visitor, [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)])
        @test visitor.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)]
        gate_visitor = Quasar.QasmGateDefVisitor(visitor, Dict{String, FreeParameter}(), Dict{String, Quasar.Qubit}(), Dict{String, Vector{Int}}(), 0, BraketSimulator.Instruction[])
        push!(gate_visitor, [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)])
        @test gate_visitor.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)]
        function_visitor = Quasar.QasmFunctionVisitor(visitor, Quasar.QasmExpression[], Quasar.QasmExpression[])
        push!(function_visitor, [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)])
        @test function_visitor.instructions == [BraketSimulator.Instruction(BraketSimulator.X(), 0), BraketSimulator.Instruction(BraketSimulator.Y(), 1)]
    end
    @testset "Sized types" begin
        for (t_string, T) in (("SizedBitVector", Quasar.SizedBitVector),
                              ("SizedInt", Quasar.SizedInt), 
                              ("SizedUInt", Quasar.SizedUInt),
                              ("SizedFloat", Quasar.SizedFloat),
                              ("SizedAngle", Quasar.SizedAngle),
                              ("SizedComplex", Quasar.SizedComplex))
            object = T(Quasar.QasmExpression(:integer_literal, 4))
            @test sprint(show, object) == t_string * "{4}"
            array_object = Quasar.SizedArray(object, (10,))
            @test sprint(show, array_object) == "SizedArray{" * t_string * "{4}, 1}"
            @test_throws BoundsError size(array_object, 1)
            @test size(array_object, 0) == 10
            new_object = T(object)
            @test sprint(show, new_object) == t_string * "{4}"
        end
        bitvector = Quasar.SizedBitVector(Quasar.QasmExpression(:integer_literal, 10))
        @test length(bitvector) == Quasar.QasmExpression(:integer_literal, 10)
        @test size(bitvector) == (Quasar.QasmExpression(:integer_literal, 10),)
    end
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
        parsed = parse_qasm(sv_adder_qasm)
        visitor = QasmProgramVisitor(Dict("a_in"=>3, "b_in"=>7))
        visitor(parsed)
        @test visitor.qubit_count == 10
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
        for (ix, c_ix) in zip(visitor.instructions, correct_instructions)
            @test ix == c_ix
        end
        @test visitor.results == BraketSimulator.Result[BraketSimulator.Probability(BraketSimulator.QubitSet(9, 5, 6, 7, 8)),
                                                        BraketSimulator.Probability(BraketSimulator.QubitSet(9)),
                                                        BraketSimulator.Probability(BraketSimulator.QubitSet(5, 6, 7, 8))]
    end
    @testset "Randomized Benchmarking" begin
        qasm = """
        qubit[2] q;
        bit[2] c;

        h q[0];
        cz q[0], q[1];
        s q[0];
        cz q[0], q[1];
        s q[0];
        z q[0];
        h q[0];
        measure q -> c;
        """
        circuit    = BraketSimulator.Circuit(qasm)
        @test circuit.instructions == [BraketSimulator.Instruction(BraketSimulator.H(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.CZ(), [0, 1]),
                                       BraketSimulator.Instruction(BraketSimulator.S(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.CZ(), [0, 1]),
                                       BraketSimulator.Instruction(BraketSimulator.S(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.Z(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.H(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.Measure(), 0),
                                       BraketSimulator.Instruction(BraketSimulator.Measure(), 1)]
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
        circuit    = BraketSimulator.Circuit(qasm) 
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        sv = 1/√2 * [-1; 0; 0; 1]
        @test simulation.state_vector ≈ sv 
        @test circuit.result_types == [BraketSimulator.Amplitude(["00", "01", "10", "11"])]
    end
    @testset "Numbers" begin
        qasm = """
        float[32] a = 1.24e-3;
        complex[float] b = 1-0.23im;
        const bit c = "0";
        bool d = false;
        complex[float] e = -0.23+2im;
        uint f = 0x123456789abcdef;
        int g = 0o010;
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["a"].val == 1.24e-3
        @test visitor.classical_defs["b"].val == 1-0.23im
        @test visitor.classical_defs["c"].val == falses(1)
        @test visitor.classical_defs["d"].val == false
        @test visitor.classical_defs["e"].val == -0.23 + 2im
        @test visitor.classical_defs["f"].val == 0x123456789abcdef
        @test visitor.classical_defs["g"].val == 0o010
    end
    @testset "Physical qubits" begin
        qasm = """
        h \$0;
        cnot \$0, \$1;
        """
        @test BraketSimulator.Circuit(qasm).instructions == [BraketSimulator.Instruction(BraketSimulator.H(), 0), BraketSimulator.Instruction(BraketSimulator.CNot(), [0, 1])]
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
        inputs = Dict("theta"=>0.2)
        parsed_circ = BraketSimulator.Circuit(qasm_str, inputs)
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
    @testset "For Loop" begin
        qasm = """
        int[8] x = 0;
        int[8] y = -100;
        int[8] ten = 10;

        for uint[8] i in [0:2:ten - 3] {
            x += i;
        }

        for int[8] i in {2, 4, 6} {
            y += i;
        }
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["x"].val == sum((0, 2, 4, 6))
        @test visitor.classical_defs["y"].val == sum((-100, 2, 4, 6))
        # without scoping { } 
        qasm = """
        int[8] x = 0;
        int[8] y = -100;
        int[8] ten = 10;

        for uint[8] i in [0:2:ten - 3] x += i;
        for int[8] i in {2, 4, 6} y += i;
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["x"].val == sum((0, 2, 4, 6))
        @test visitor.classical_defs["y"].val == sum((-100, 2, 4, 6))
    end
    @testset "While Loop" begin
        qasm = """
        int[8] x = 0;
        int[8] i = 0;

        while (i < 7) {
            x += i;
            i += 1;
        }
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["x"].val == sum(0:6)
        # without scoping { }
        qasm = """
        int[8] i = 0;

        while (i < 7) i += 1;
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["i"].val == 7
    end
    @testset "Break and continue" begin
        @testset for str in ("continue;", "{ continue; }")
            qasm = """
            int[8] x = 0;
            for int[8] i in [0: 3] {
                if (mod(i, 2) == 0) $str
                x += i;
            }
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            visitor(parsed)
            @test visitor.classical_defs["x"].val == 4
            qasm = """
            int[8] x = 0;
            int[8] i = 0;

            while (i < 7) {
                i += 1;
                if (i == 5) $str
                x += i;
            }
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            visitor(parsed)
            @test visitor.classical_defs["x"].val == sum(0:7) - 5
        end
        @testset for str in ("break;", "{ break; }")
            qasm = """
            int[8] x = 0;
            for int[8] i in [0: 3] {
                x += i;
                if (mod(i, 2) == 1) $str
            }
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            visitor(parsed)
            @test visitor.classical_defs["x"].val == 1
            qasm = """
            int[8] x = 0;
            int[8] i = 0;

            while (i < 7) {
                x += i;
                i += 1;
                if (i == 5) $str
            }
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            visitor(parsed)
            @test visitor.classical_defs["x"].val == sum(0:4)
        end
        @testset for str in ("continue", "break")
            qasm = """
            int[8] x = 0;
            int[8] i = 0;
            $str;
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            @test_throws Quasar.QasmVisitorError("$str statement encountered outside a loop.") visitor(parsed)
        end
    end
    @testset "Builtin functions" begin
        qasm = """
            const float[64] arccos_result = arccos(1);
            const float[64] arcsin_result = arcsin(1);
            const float[64] arctan_result = arctan(1);
            const int[64] ceiling_result = ceiling(π);
            const float[64] cos_result = cos(1);
            const float[64] exp_result = exp(2);
            const int[64] floor_result = floor(π);
            const float[64] log_result = log(ℇ);
            const int[64] mod_int_result = mod(4, 3);
            const float[64] mod_float_result = mod(5.2, 2.5);
            const int[64] popcount_bool_result = popcount(true);
            const int[64] popcount_bit_result = popcount("1001110");
            const int[64] popcount_int_result = popcount(78);
            const int[64] popcount_uint_result = popcount(0x78);
            bit[2] bitvec;
            bitvec = "11";
            const int[64] popcount_bitvec_result = popcount(bitvec);
            // parser gets confused by pow
            // const int[64] pow_int_result = pow(3, 3);
            // const float[64] pow_float_result = pow(2.5, 2.5);
            // add rotl, rotr
            const float[64] sin_result = sin(1);
            const float[64] sqrt_result = sqrt(2);
            const float[64] tan_result = tan(1);
            """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["arccos_result"].val == acos(1)
        @test visitor.classical_defs["arcsin_result"].val == asin(1)
        @test visitor.classical_defs["arctan_result"].val == atan(1)
        @test visitor.classical_defs["ceiling_result"].val == 4
        @test visitor.classical_defs["cos_result"].val == cos(1)
        @test visitor.classical_defs["exp_result"].val == exp(2)
        @test visitor.classical_defs["floor_result"].val == 3
        @test visitor.classical_defs["log_result"].val == 1.0
        @test visitor.classical_defs["mod_int_result"].val == 1
        @test visitor.classical_defs["mod_float_result"].val == mod(5.2, 2.5)
        @test visitor.classical_defs["popcount_bool_result"].val == 1
        @test visitor.classical_defs["popcount_bit_result"].val == 4
        @test visitor.classical_defs["popcount_int_result"].val == 4
        @test visitor.classical_defs["popcount_uint_result"].val == 4
        @test visitor.classical_defs["popcount_bitvec_result"].val == 2
        @test visitor.classical_defs["sin_result"].val == sin(1)
        @test visitor.classical_defs["sqrt_result"].val == sqrt(2)
        @test visitor.classical_defs["tan_result"].val == tan(1)

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
            parsed_circ = BraketSimulator.Circuit(qasm, inputs) 
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
    @testset "Reset" begin
        qasm = """
        qubit[4] q;
        x q[0];
        reset q[0];
        """
        @test_warn "reset expression encountered -- currently `reset` is a no-op" parse_qasm(qasm)
    end
    @testset "Gate call missing/extra args" begin
        qasm = """
        qubit[4] q;
        rx q[0];
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        @test_throws Quasar.QasmVisitorError("gate rx requires arguments but none were provided.") visitor(parsed)
        qasm = """
        qubit[4] q;
        x(0.1) q[0];
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        @test_throws Quasar.QasmVisitorError("gate x does not accept arguments but arguments were provided.") visitor(parsed)
    end
    #=@testset "Adjoint Gradient pragma" begin
        qasm = """
        input float theta;
        qubit[4] q;
        rx(theta) q[0];
        #pragma braket result adjoint_gradient expectation(-6 * y(q[0]) @ i(q[1]) + 0.75 * y(q[2]) @ z(q[3])) theta
        """
        circuit = BraketSimulator.Circuit(qasm, Dict("theta"=>0.1))
        θ       = FreeParameter("theta")
        obs     = -2 * BraketSimulator.Observables.Y() * (3 * BraketSimulator.Observables.I()) + 0.75 * BraketSimulator.Observables.Y() * BraketSimulator.Observables.Z()
        @test circuit.result_types[1].observable == obs
        @test circuit.result_types[1].targets == [QubitSet([0, 1]), QubitSet([2, 3])]
        @test circuit.result_types[1].parameters == ["theta"]
    end=#
    @testset "Assignment operators" begin
        qasm = """
        int[16] x;
        bit[4] xs;

        x = 0;
        xs = "0000";

        x += 1; // 1
        x *= 2; // 2
        x /= 2; // 1
        x -= 5; // -4

        xs[2:] |= "11";
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["x"].val  == -4
        @test visitor.classical_defs["xs"].val == BitVector([0,0,1,1])
    end
    @testset "Bit operators" begin
        qasm = """
        bit[4] and;
        bit[4] or;
        bit[4] xor;
        bit[4] lshift;
        bit[4] rshift;
        bit[4] flip;
        bit gt;
        bit lt;
        bit ge;
        bit le;
        bit eq;
        bit neq;
        bit not;
        bit not_zero;

        bit[4] x = "0101";
        bit[4] y = "1100";

        and = x & y;
        or  = x | y;
        xor = x ^ y;
        lshift = x << 2;
        rshift = y >> 2;
        flip = ~x;
        gt = x > y;
        lt = x < y;
        ge = x >= y;
        le = x <= y;
        eq = x == y;
        neq = x != y;
        not = !x;
        not_zero = !(x << 4);
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["x"].val == BitVector([0,1,0,1])
        @test visitor.classical_defs["y"].val == BitVector([1,1,0,0])
        @test visitor.classical_defs["and"].val == BitVector([0,1,0,0])
        @test visitor.classical_defs["or"].val == BitVector([1,1,0,1])
        @test visitor.classical_defs["xor"].val == BitVector([1,0,0,1])
        @test visitor.classical_defs["lshift"].val == BitVector([0,1,0,0])
        @test visitor.classical_defs["rshift"].val == BitVector([0,0,1,1])
        @test visitor.classical_defs["flip"].val == BitVector([1,0,1,0])
        @test visitor.classical_defs["gt"].val == false
        @test visitor.classical_defs["lt"].val == true
        @test visitor.classical_defs["ge"].val == false
        @test visitor.classical_defs["le"].val == true
        @test visitor.classical_defs["eq"].val == false
        @test visitor.classical_defs["neq"].val == true
        @test visitor.classical_defs["not"].val == false
        @test visitor.classical_defs["not_zero"].val == true
    end
    @testset "End statement" begin
        qasm = """
        z \$0;
        x \$1;
        h \$2;
        y \$3;
        end;
        y \$2;
        h \$3;
        """
        circ = BraketSimulator.Circuit(qasm)
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0),
                                    BraketSimulator.Instruction(BraketSimulator.X(), 1),
                                    BraketSimulator.Instruction(BraketSimulator.H(), 2),
                                    BraketSimulator.Instruction(BraketSimulator.Y(), 3)]
    end
    @testset "Switch/case" begin
        qasm = """
        input int[8] x;
        switch (x + 1) {
            case 0b00 {}
            default { z \$0; }
        }
        """
        circ = BraketSimulator.Circuit(qasm, Dict("x"=> -1))
        @test isempty(circ.instructions)
        circ = BraketSimulator.Circuit(qasm, Dict("x"=> 0))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        qasm = """
        input int[8] x;
        switch (x) { case 0 {} case 1, 2 { z \$0; }  }
        """
        circ = BraketSimulator.Circuit(qasm, Dict("x"=> 0))
        @test isempty(circ.instructions)
        circ = BraketSimulator.Circuit(qasm, Dict("x"=> 1))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        circ = BraketSimulator.Circuit(qasm, Dict("x"=> 2))
        @test circ.instructions == [BraketSimulator.Instruction(BraketSimulator.Z(), 0)]
        qasm = """
        input int[8] x;
        switch (x) { case 0 {} default { z \$0; } default { x \$0; } }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.Circuit(qasm, Dict("x"=>0))
        qasm = """
        input int[8] x;
        switch (x) { default { z \$0; } case 0 {} }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.Circuit(qasm, Dict("x"=>0))
        qasm = """
        input int[8] x;
        switch (x) { case 0 { z \$0; } true {} }
        """
        @test_throws Quasar.QasmParseError BraketSimulator.Circuit(qasm, Dict("x"=>0))
    end
    @testset "If/Else" begin
        qasm = """
        int[8] two = 2;
        bit[3] m = "000";

        if (two + 1) {
            m[0] = 1;
        } else {
            m[1] = 1;
        }

        if (!bool(two - 2)) {
            m[2] = 1;
        }
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        bit_vec = BitVector([1,0,1])
        @test visitor.classical_defs["m"].val == bit_vec
        qasm = """
        int[8] two = -2;
        bit[3] m = "000";

        if (two + 1) {
            m[0] = 1;
        } else {
            m[1] = 1;
        }

        if (!bool(two - 2)) {
            m[2] = 1;
        }
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        bit_vec = BitVector([0,1,1])
        @test visitor.classical_defs["m"].val == bit_vec
        qasm = """
        input bool which_gate;
        qubit q;

        if (which_gate) {
            x q;
        } else {
            y q;
        }
        """
        for (flag, which_gate) in ((true, BraketSimulator.X()), (false, BraketSimulator.Y()))
            circuit = BraketSimulator.Circuit(qasm, Dict("which_gate"=>flag))
            @test circuit.instructions == [BraketSimulator.Instruction(which_gate, 0)]
        end
    end
    @testset "Global gate control" begin
        qasm = """
        qubit q1;
        qubit q2;

        h q1;
        h q2;
        ctrl @ s q1, q2;
        """
        circuit    = BraketSimulator.Circuit(qasm) 
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ [0.5, 0.5, 0.5, 0.5im]
    end
    @testset "Simple Pow" begin
        custom_qasm = """
        qubit q1;
        qubit q2;
        h q1;
        h q2;
        
        pow(1/2) @ z q1; // s
        """
        standard_qasm = """
        qubit q1;
        qubit q2;
        h q1;
        h q2;
        
        s q1; // s
        """
        canonical_ixs = [BraketSimulator.Instruction(BraketSimulator.H(), 0), BraketSimulator.Instruction(BraketSimulator.H(), 1), BraketSimulator.Instruction(BraketSimulator.S(), 0)]
        canonical_simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(canonical_simulation, canonical_ixs)
        @testset "$title" for (title, qasm) in (("standard", standard_qasm), ("custom", custom_qasm))
            circuit    = BraketSimulator.Circuit(qasm)
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
            circuit    = BraketSimulator.Circuit(qasm) 
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
        circuit    = BraketSimulator.Circuit(qasm) 
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
        circuit    = BraketSimulator.Circuit(qasm)
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
        circuit   = BraketSimulator.Circuit(qasm) 
        collapsed = prod(BraketSimulator.matrix_rep(ix.operator) for ix in circuit.instructions)
        @test collapsed ≈ diagm(ones(ComplexF64, 2^BraketSimulator.qubit_count(circuit)))
    end
    @testset "Noise" begin
        qasm = """
        qubit[2] qs;

        #pragma braket noise bit_flip(.5) qs[1]
        #pragma braket noise phase_flip(.5) qs[0]
        #pragma braket noise pauli_channel(.1, .2, .3) qs[0]
        #pragma braket noise depolarizing(.5) qs[0]
        #pragma braket noise two_qubit_depolarizing(.9) qs
        #pragma braket noise two_qubit_depolarizing(.7) qs[1], qs[0]
        #pragma braket noise two_qubit_dephasing(.6) qs
        #pragma braket noise amplitude_damping(.2) qs[0]
        #pragma braket noise generalized_amplitude_damping(.2, .3)  qs[1]
        #pragma braket noise phase_damping(.4) qs[0]
        #pragma braket noise kraus([[0.9486833im, 0], [0, 0.9486833im]], [[0, 0.31622777], [0.31622777, 0]]) qs[0]
        #pragma braket noise kraus([[0.9486832980505138, 0, 0, 0], [0, 0.9486832980505138, 0, 0], [0, 0, 0.9486832980505138, 0], [0, 0, 0, 0.9486832980505138]], [[0, 0.31622776601683794, 0, 0], [0.31622776601683794, 0, 0, 0], [0, 0, 0, 0.31622776601683794], [0, 0, 0.31622776601683794, 0]]) qs[{1, 0}]
        """
        circuit   = BraketSimulator.Circuit(qasm) 
        inst_list = [BraketSimulator.Instruction(BraketSimulator.BitFlip(0.5), [1]),
                     BraketSimulator.Instruction(BraketSimulator.PhaseFlip(0.5), [0]),
                     BraketSimulator.Instruction(BraketSimulator.PauliChannel(0.1, 0.2, 0.3), [0]),
                     BraketSimulator.Instruction(BraketSimulator.Depolarizing(0.5), [0]),
                     BraketSimulator.Instruction(BraketSimulator.TwoQubitDepolarizing(0.9), [0, 1]),
                     BraketSimulator.Instruction(BraketSimulator.TwoQubitDepolarizing(0.7), [1, 0]),
                     BraketSimulator.Instruction(BraketSimulator.TwoQubitDephasing(0.6), [0, 1]),
                     BraketSimulator.Instruction(BraketSimulator.AmplitudeDamping(0.2), [0]),
                     BraketSimulator.Instruction(BraketSimulator.GeneralizedAmplitudeDamping(0.2, 0.3), [1]),
                     BraketSimulator.Instruction(BraketSimulator.PhaseDamping(0.4), [0]),
                     BraketSimulator.Instruction(BraketSimulator.Kraus([[0.9486833im 0; 0 0.9486833im], [0 0.31622777; 0.31622777 0]]), [0]),
                     BraketSimulator.Instruction(BraketSimulator.Kraus([diagm(fill(√0.9, 4)), √0.1*kron([1.0 0.0; 0.0 1.0], [0.0 1.0; 1.0 0.0])]), [1, 0]),
                    ] 
        @testset "Operator $(typeof(ix.operator)), target $(ix.target)" for (cix, ix) in zip(circuit.instructions, inst_list)
            @test cix.operator == ix.operator
            @test cix.target == ix.target
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
            circuit = BraketSimulator.Circuit(qasm)
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
            circuit = BraketSimulator.Circuit(qasm) 
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
            circuit = BraketSimulator.Circuit(qasm) 
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
        qasm = """
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
        circuit    = BraketSimulator.Circuit(qasm)
        simulation = BraketSimulator.StateVectorSimulator(3, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ [0, 0, 0, 0, 0.70710678, 0, 0, 0.70710678]
    end
    @testset "Output" begin
        qasm = """
               output int[8] out_int;
               """
        parsed  = parse_qasm(qasm)
        visitor = Quasar.QasmProgramVisitor()
        @test_throws Quasar.QasmVisitorError("Output not supported.") visitor(parsed)
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
        circuit    = BraketSimulator.Circuit(qasm)
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
        circuit        = BraketSimulator.Circuit(with_verbatim)
        sim_w_verbatim = BraketSimulator.StateVectorSimulator(2, 0)
        pop!(circuit.instructions) 
        pop!(circuit.instructions) 
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
        circuit         = BraketSimulator.Circuit(without_verbatim)
        sim_wo_verbatim = BraketSimulator.StateVectorSimulator(2, 0)
        pop!(circuit.instructions) 
        pop!(circuit.instructions) 
        BraketSimulator.evolve!(sim_wo_verbatim, circuit.instructions)

        @test sim_w_verbatim.state_vector ≈ sim_wo_verbatim.state_vector
    end
    @testset "Subroutine" begin
        qasm = """
        const int[8] n = 4;
        def parity(bit[n] cin) -> bit {
          bit c = false;
          for int[8] i in [0: n - 1] {
            c ^= cin[i];
          }
          return c;
        }

        bit[n] c = "1011";
        bit p = parity(c);
        """
        parsed  = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["p"].val == true
    end
    @testset "Void subroutine" begin
        qasm = """
               def flip(qubit q) {
                 x q;
               }
               qubit[2] qs;
               flip(qs[0]);
               """
        circuit    = BraketSimulator.Circuit(qasm)
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        BraketSimulator.evolve!(simulation, circuit.instructions)
        @test simulation.state_vector ≈ [0, 0, 1, 0]
    end
    @testset "Array ref subroutine" begin
        qasm = """
        int[16] total_1;
        int[16] total_2;

        def sum(readonly array[int[8], #dim = 1] arr) -> int[16] {
            int[16] size = sizeof(arr);
            int[16] x = 0;
            for int[8] i in [0:size - 1] {
                x += arr[i];
            }
            return x;
        }

        array[int[8], 5] array_1 = {1, 2, 3, 4, 5};
        array[int[8], 10] array_2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        total_1 = sum(array_1);
        total_2 = sum(array_2);
        """
        parsed = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["total_1"].val == 15
        @test visitor.classical_defs["total_2"].val == 55
        @test visitor.classical_defs["array_1"].val == [1, 2, 3, 4, 5]
        @test visitor.classical_defs["array_2"].val == collect(1:10)
    end
    @testset "Const array argument" begin
        qasm_string = """
        qubit q;

        def sum(const array[int[8], #dim = 1] arr) -> int {
            int size = sizeof(arr);
            int x = 0;
            for int i in [0:size - 1] {
                x += arr[i];
            }
            return x;
        }

        array[int, 10] arr = {9, 3, 6, 2, 2, 4, 3, 1, 12, 7};
        int s = sum(arr);
        rx(s*π/4) q;
        """
        parsed = parse_qasm(qasm_string)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["s"].val == sum([ 9, 3, 6, 2, 2, 4, 3, 1, 12, 7 ])
    end
    @testset "Array ref subroutine with mutation" begin
        qasm = """
        def mutate_array(mutable array[int[8], #dim = 1] arr) {
            int[16] size = sizeof(arr);
            for int[8] i in [0:size - 1] {
                arr[i] = 0;
            }
        }

        array[int[8], 5] array_1 = {1, 2, 3, 4, 5};
        array[int[8], 10] array_2 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        array[int[8], 10] array_3 = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        mutate_array(array_1);
        mutate_array(array_2);
        mutate_array(array_3[4:2:-1]);
        """
        parsed = parse_qasm(qasm)
        visitor = QasmProgramVisitor()
        visitor(parsed)
        @test visitor.classical_defs["array_1"].val == zeros(Int, 5) 
        @test visitor.classical_defs["array_2"].val == zeros(Int, 10) 
        @test visitor.classical_defs["array_3"].val == [1, 2, 3, 4, 0, 6, 0, 8, 0, 10]
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
            circuit    = BraketSimulator.Circuit(qasm) 
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
        circuit    = BraketSimulator.Circuit(qasm)
        program    = BraketSimulator.Program(circuit)
        simulation = BraketSimulator.StateVectorSimulator(2, 0)
        @test_throws BraketSimulator.ValidationError BraketSimulator.simulate(simulation, program, 2, 0)
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
    @testset "Qubits with variable as size" begin
        qasm_string = """
        OPENQASM 3.0;

        def ghz(int[32] n) {
            h q[0];
            for int i in [0:n - 1] {
                cnot q[i], q[i + 1];
            }
        }

        int[32] n = 5;
        bit[n + 1] c;
        qubit[n + 1] q;

        ghz(n);

        c = measure q;
        """
        parsed = parse_qasm(qasm_string)
        visitor = Quasar.QasmProgramVisitor()
        visitor(parsed)
        @test Quasar.qubit_count(visitor) == 6
    end
    @testset "String inputs" begin
        qasm_string = """
        const int[8] n = 4;
        input bit[n] x;

        qubit q;

        def parity(bit[n] cin) -> bit {
          bit c = false;
          for int[8] i in [0: n - 1] {
            c ^= cin[i];
          }
          return c;
        }

        if(parity(x)) {
            x q;
        } else {
            i q;
        }
        """
        circ = BraketSimulator.Circuit(qasm_string, Dict("x"=>"1011"))
    end
    @testset "GRCS 16" begin
        simulator = BraketSimulator.StateVectorSimulator(16, 0)
        circuit   = BraketSimulator.Circuit(joinpath(@__DIR__, "grcs_16.qasm")) 
        BraketSimulator.evolve!(simulator, circuit.instructions)
        probs     = BraketSimulator.probabilities(simulator)        
        @test first(probs) ≈ 0.0000062 atol=1e-7
    end
    @testset "Include" begin
        mktempdir() do dir_path
            inc_path = joinpath(dir_path, "new_gate.inc")
            write(inc_path, """\ngate u3(θ, ϕ, λ) q { gphase(-(ϕ+λ)/2); U(θ, ϕ, λ) q; }\n""")
            qasm = """
            OPENQASM 3;
            include "$inc_path";
            """
            parsed  = parse_qasm(qasm)
            visitor = QasmProgramVisitor()
            visitor(parsed)
            @test haskey(visitor.gate_defs, "u3")
        end
    end
end
