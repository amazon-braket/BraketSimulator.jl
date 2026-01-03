using Test,
      Logging,
      Statistics,
      LinearAlgebra,
      BraketSimulator

using Braket
using Braket: I, name

@testset "Basic integration of local simulators with BraketSimulator.jl" begin
    @testset "Simulator $sim_type" for (sim_type, rt) in (
        ("braket_sv_v2", Braket.StateVector),
        ("braket_dm_v2", Braket.DensityMatrix),
    )
        d = LocalSimulator(sim_type)
        @test d.backend == sim_type
        c = Braket.Circuit()
        Braket.H(c, 0, 1, 2)
        Braket.Rx(c, 0, 1, 2, 0.5)
        rt(c)
        if sim_type == "braket_sv_v2"
            Braket.Amplitude(c, ["000", "111"])
        end
        with_logger(NullLogger()) do
            r = d(c, shots = 0)
        end
    end
end

@testset "BraketSimulator to Braket Conversion" begin
@testset "Simple Circuit" begin
  # Create a new circuit
  circuit = Braket.Circuit()

  # Add instructions to the circuit
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 0))          # Apply Hadamard to qubit 0
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.CZ(), 0, 1))      # Apply controlled-Z between qubit 0 and qubit 1
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.Measure(), 0))    # Measure qubit 0
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.Measure(), 1))    # Measure qubit 1

  # Convert the circuit to QASM format
  qasm = """
       qubit[2] q;
       bit[2] c;

       h q[0];
       cz q[0], q[1];
       measure q -> c;
       """

  # Convert the QASM representation back to a Braket circuit
  converted_circuit = convert(Braket.Circuit, BraketSimulator.Circuit(qasm))

  # Assert that the original circuit and the converted circuit are equivalent
  @test circuit == converted_circuit
end

@testset "Active Reset and Mid-Circuit Measurement Circuit" begin
  # Blocked by Braket.jl/issues/99 and PR #46
  # BLocked by Reset Instruction missing in Braket
  # Create a new circuit using Braket
  circuit = Braket.Circuit()

  # Add instructions to the circuit
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.X(), 2))          # Apply X gate to qubit 2
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 2))          # Apply Hadamard gate to qubit 2
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 1))          # Apply Hadamard gate to qubit 1
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.CNot(), 1, 2))    # Apply CNot gate with qubit 1 as control and qubit 2 as target
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 1))          # Apply Hadamard gate to qubit 1
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.Measure(), 1))    # Measure qubit 1 and store result in bit 0
  # Uncomment when Reset is added to Braket
  # Braket.add_instruction!(circuit, Braket.Instruction(Braket.Reset(), 1))      # Reset qubit 1
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 1))          # Apply Hadamard gate to qubit 1 again
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.CNot(), 1, 2))    # Apply CNot gate again
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.H(), 1))          # Apply Hadamard gate to qubit 1 again
  Braket.add_instruction!(circuit, Braket.Instruction(Braket.Measure(), 1))    # Measure qubit 1 again and store result in bit 0
  # Uncomment when Reset is added to Braket
  # Braket.add_instruction!(circuit, Braket.Instruction(Braket.Reset(), 1))      # Reset qubit 1 again

  # Convert the circuit to QASM format
  qasm = """
  OPENQASM 3.0;
  bit[2] b;
  qubit[3] q;
  x q[2];
  h q[2];
  h q[1];
  cnot q[1], q[2];
  h q[1];
  b[0] = measure q[1];
  // reset q[1]; // TODO add when Braket supports reset
  h q[1];
  cnot q[1], q[2];
  h q[1];
  b[0] = measure q[1];
  // reset q[1]; // TODO add when Braket supports reset
  """

  # Convert the QASM representation back to a Braket circuit
  converted_circuit = convert(Braket.Circuit, BraketSimulator.Circuit(qasm))

  # Assert that the original circuit and the converted circuit are equivalent
  @test circuit == converted_circuit
end
end

@testset "Type conversions" begin
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.Measure(2))) == Braket.Measure(2)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.Reset())) == Braket.Reset()
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.Barrier())) == Braket.Barrier()
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.Delay(Braket.Nanosecond(2)))) == Braket.Delay(Braket.Nanosecond(2))
    angle1 = 0.2
    angle2 = 0.1
    angle3 = π
    prob = 0.015
    prob2 = 0.1
    prob3 = 0.002
    gamma = 0.23
    prob_dict = Dict("XX"=>prob, "YY"=>prob2, "ZZ"=>prob3)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.MS(angle1, angle2, angle3))) == Braket.MS(angle1, angle2, angle3)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.U(angle1, angle2, angle3))) == Braket.U(angle1, angle2, angle3)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.PRx(angle1, angle2))) == Braket.PRx(angle1, angle2)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.GeneralizedAmplitudeDamping(prob, gamma))) == Braket.GeneralizedAmplitudeDamping(prob, gamma)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.PauliChannel(prob, prob3, prob3))) == Braket.PauliChannel(prob, prob3, prob3)
    @test convert(Braket.Operator, convert(BraketSimulator.Operator, Braket.MultiQubitPauliChannel{2}(prob_dict))) == Braket.MultiQubitPauliChannel{2}(prob_dict)
    @test convert(Braket.AbstractProgramResult, convert(BraketSimulator.AbstractProgramResult, Braket.IR.StateVector("statevector"))) == Braket.IR.StateVector("statevector")
    @test qubit_count(BraketSimulator.Observables.X()) == 1
    @test qubit_count(BraketSimulator.Observables.Y()) == 1
    @test qubit_count(BraketSimulator.Observables.Z()) == 1
    @test qubit_count(BraketSimulator.Observables.TensorProduct(["x", "h", "y"])) == 3
end

@testset "Correctness" begin
    PURE_DEVICE = LocalSimulator("braket_sv_v2")
    NOISE_DEVICE = LocalSimulator("braket_dm_v2")
    ALL_DEVICES = [PURE_DEVICE, NOISE_DEVICE]
    PURE_DEVICES = [PURE_DEVICE]
    NOISE_DEVICES = [NOISE_DEVICE]
    SHOT_LIST = (0, 8000)

    # looser tolerance bounds here to account
    # for differences in `np.allclose` vs `isapprox`
    get_tol(shots::Int) = return (
        shots > 0 ? Dict("atol" => 0.2, "rtol" => 0.25) : Dict("atol" => 0.01, "rtol" => 0)
    )

    bell_circ() = Braket.Circuit([(Braket.H, 0), (Braket.CNot, 0, 1)])
    three_qubit_circuit(
        θ::Float64,
        ϕ::Float64,
        φ::Float64,
        obs::Braket.Observables.Observable,
        obs_targets::Vector{Int},
    ) = Braket.Circuit([
        (Braket.Rx, 0, θ),
        (Braket.Rx, 1, ϕ),
        (Braket.Rx, 2, φ),
        (Braket.CNot, 0, 1),
        (Braket.CNot, 1, 2),
        (Braket.Variance, obs, obs_targets),
        (Braket.Expectation, obs, obs_targets),
    ])

    @inline function variance_expectation_sample_result(
        res::Braket.GateModelQuantumTaskResult,
        shots::Int,
        expected_var::Float64,
        expected_mean::Float64,
        expected_eigs::Vector{Float64},
    )
        tol = get_tol(shots)
        variance = res.values[1]
        expectation = res.values[2]
        if shots > 0
            samples = res.values[3]
            sign_fix(x) = (iszero(x) || abs(x) < 1e-12) ? 0.0 : x
            fixed_samples = sort(collect(unique(sign_fix, samples)))
            fixed_eigs = sort(collect(unique(sign_fix, expected_eigs)))
            @test isapprox(
                sort(fixed_samples),
                sort(fixed_eigs),
                rtol = tol["rtol"],
                atol = tol["atol"],
            )
            @test isapprox(
                mean(samples),
                expected_mean,
                rtol = tol["rtol"],
                atol = tol["atol"],
            )
            @test isapprox(
                var(samples),
                expected_var,
                rtol = tol["rtol"],
                atol = tol["atol"],
            )
        end
        @test isapprox(expectation, expected_mean, rtol = tol["rtol"], atol = tol["atol"])
        @test isapprox(variance, expected_var, rtol = tol["rtol"], atol = tol["atol"])
    end

    @testset "Local Braket Simulator" begin
        @testset for (backend, device_name) in [
            ("braket_sv_v2", "StateVectorSimulator"),
            ("braket_dm_v2", "DensityMatrixSimulator"),
        ]
            local_simulator_device = LocalSimulator(backend)
            @test Braket.name(local_simulator_device) == device_name
        end
        @testset "Device $DEVICE, shots $SHOTS" for DEVICE in ALL_DEVICES,
            SHOTS in SHOT_LIST
            with_logger(NullLogger()) do
                if SHOTS > 0
                    @testset "qubit ordering" begin
                        device = DEVICE
                        state_110 = Braket.Circuit([(Braket.X, 0), (Braket.X, 1), (Braket.I, 2)])
                        state_001 = Braket.Circuit([(Braket.I, 0), (Braket.I, 1), (Braket.X, 2)])
                        @testset for (state, most_com) in
                                     ((state_110, "110"), (state_001, "001"))
                            tasks = (state, ir(state, Val(:JAQCD)), ir(state, Val(:OpenQASM)))
                            @testset for task in tasks
                                res = result(device(task, shots = SHOTS))
                                mc  = argmax(res.measurement_counts)
                                @test mc == most_com
                            end
                        end
                    end

                    @testset "Bell pair nonzero shots" begin
                        circuit = bell_circ()
                        circuit(Braket.Expectation, Braket.Observables.H() * Braket.Observables.X(), [0, 1])
                        circuit(Braket.Sample, Braket.Observables.H() * Braket.Observables.X(), [0, 1])
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        @testset for task in tasks
                            device = DEVICE
                            res = result(device(task; shots = SHOTS))
                            @test length(res.result_types) == 2
                            @test 0.6 <
                                  res[Braket.Expectation(Braket.Observables.H() * Braket.Observables.X(), [0, 1])] <
                                  0.8
                            @test length(
                                res[Braket.Sample(Braket.Observables.H() * Braket.Observables.X(), [0, 1])],
                            ) == SHOTS
                        end
                    end
                end
                @testset "Bell pair full probability" begin
                    circuit = bell_circ()
                    circuit(Braket.Probability)
                    tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    tol = get_tol(SHOTS)
                    @testset for task in tasks
                        device = DEVICE
                        res = result(device(task, shots = SHOTS))
                        @test length(res.result_types) == 1
                        @test isapprox(
                            res[Probability()],
                            [0.5, 0.0, 0.0, 0.5],
                            rtol = tol["rtol"],
                            atol = tol["atol"],
                        )
                    end
                end
                @testset "Bell pair marginal probability" begin
                    circuit = bell_circ()
                    circuit(Braket.Probability, 0)
                    tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    tol = get_tol(SHOTS)
                    @testset for task in tasks
                        device = DEVICE
                        res = result(device(task, shots = SHOTS))
                        @test length(res.result_types) == 1
                        @test isapprox(
                            res[Braket.Probability(0)],
                            [0.5, 0.5],
                            rtol = tol["rtol"],
                            atol = tol["atol"],
                        )
                    end
                end
                @testset "Result types x x y" begin
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    obs_targets = [0, 2]
                    expected_mean = sin(θ) * sin(ϕ) * sin(φ)
                    expected_var =
                        (
                            8 * sin(θ)^2 * cos(2φ) * sin(ϕ)^2 - cos(2(θ - ϕ)) - cos(2(θ + ϕ)) +
                            2 * cos(2θ) +
                            2 * cos(2ϕ) +
                            14
                        ) / 16
                    expected_eigs = [-1.0, 1.0]
                    device = DEVICE
                    shots = SHOTS
                    @testset "Obs $obs" for obs in (
                        Braket.Observables.X() * Braket.Observables.Y(),
                        Braket.Observables.HermitianObservable(kron([0 1; 1 0], [0 -im; im 0])),
                    )
                        circuit = three_qubit_circuit(θ, ϕ, φ, obs, obs_targets)
                        shots > 0 && circuit(Braket.Sample, obs, obs_targets)
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        for task in tasks
                            res = result(device(task, shots = shots))
                            variance_expectation_sample_result(
                                res,
                                shots,
                                expected_var,
                                expected_mean,
                                expected_eigs,
                            )
                        end
                    end
                end
                @testset "Result types z x h x y" begin
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    obs = Braket.Observables.Z() * Braket.Observables.H() * Braket.Observables.Y()
                    obs_targets = [0, 1, 2]
                    circuit = three_qubit_circuit(θ, ϕ, φ, obs, obs_targets)
                    expected_mean = -(cos(φ) * sin(ϕ) + sin(φ) * cos(θ)) / √2
                    expected_var =
                        (
                            3 + cos(2ϕ) * cos(φ)^2 - cos(2θ) * sin(φ)^2 -
                            2 * cos(θ) * sin(ϕ) * sin(2φ)
                        ) / 4
                    expected_eigs = [-1.0, 1.0]
                    device = DEVICE
                    shots = SHOTS
                    @testset for obs in (
                        Braket.Observables.Z() * Braket.Observables.H() * Braket.Observables.Y(),
                    )
                        circuit = three_qubit_circuit(θ, ϕ, φ, obs, obs_targets)
                        shots > 0 && circuit(Braket.Sample, obs, obs_targets)
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        for task in tasks
                            res = result(device(task, shots = shots))
                            variance_expectation_sample_result(
                                res,
                                shots,
                                expected_var,
                                expected_mean,
                                expected_eigs,
                            )
                        end
                    end
                end
                @testset "Result types z x z" begin
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    obs_targets = [0, 2]
                    expected_mean = 0.849694136476246
                    expected_var = 0.27801987443788634
                    expected_eigs = [-1.0, 1.0]
                    device = DEVICE
                    shots = SHOTS
                    @testset for obs in (
                        Braket.Observables.Z() * Braket.Observables.Z(),
                        Braket.Observables.HermitianObservable(kron([1 0; 0 -1], [1 0; 0 -1])),
                    )
                        circuit = three_qubit_circuit(θ, ϕ, φ, obs, obs_targets)
                        shots > 0 && circuit(Braket.Sample, obs, obs_targets)
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        @testset for task in tasks
                            res = result(device(task, shots = shots))
                            variance_expectation_sample_result(
                                res,
                                shots,
                                expected_var,
                                expected_mean,
                                expected_eigs,
                            )
                        end
                    end
                end
                @testset "($DEVICE, $SHOTS) Result types tensor {i,y,z,Hermitian} x Hermitian" begin
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    ho_mat = [
                        -6 2+im -3 -5+2im
                        2-im 0 2-im -5+4im
                        -3 2+im 0 -4+3im
                        -5-2im -5-4im -4-3im -6
                    ]
                    @test ishermitian(ho_mat)
                    ho = Braket.Observables.HermitianObservable(ComplexF64.(ho_mat))
                    ho_mat2 = [1 2; 2 4]
                    ho2 = Braket.Observables.HermitianObservable(ComplexF64.(ho_mat2))
                    ho_mat3 = [-6 2+im; 2-im 0]
                    ho3 = Braket.Observables.HermitianObservable(ComplexF64.(ho_mat3))
                    ho_mat4 = kron([1 0; 0 1], [-6 2+im; 2-im 0])
                    ho4 = Braket.Observables.HermitianObservable(ComplexF64.(ho_mat4))
                    meani = -5.7267957792059345
                    meany = 1.4499810303182408
                    meanz =
                        0.5 * (
                            -6 * cos(θ) * (cos(φ) + 1) -
                            2 * sin(φ) * (cos(θ) + sin(ϕ) - 2 * cos(ϕ)) +
                            3 * cos(φ) * sin(ϕ) +
                            sin(ϕ)
                        )
                    meanh = -4.30215023196904
                    meanii = -5.78059066879935

                    vari = 43.33800156673375
                    vary = 74.03174647518193
                    varz =
                        (
                            1057 - cos(2ϕ) + 12 * (27 + cos(2ϕ)) * cos(φ) -
                            2 * cos(2φ) * sin(ϕ) * (16 * cos(ϕ) + 21 * sin(ϕ)) + 16 * sin(2ϕ) -
                            8 * (-17 + cos(2ϕ) + 2 * sin(2ϕ)) * sin(φ) -
                            8 * cos(2θ) * (3 + 3 * cos(φ) + sin(φ))^2 -
                            24 * cos(ϕ) * (cos(ϕ) + 2 * sin(ϕ)) * sin(2φ) -
                            8 *
                            cos(θ) *
                            (
                                4 *
                                cos(ϕ) *
                                (4 + 8 * cos(φ) + cos(2φ) - (1 + 6 * cos(φ)) * sin(φ)) +
                                sin(ϕ) *
                                (15 + 8 * cos(φ) - 11 * cos(2φ) + 42 * sin(φ) + 3 * sin(2φ))
                            )
                        ) / 16
                    varh = 370.71292282796804
                    varii = 6.268315532585994

                    i_array = [1 0; 0 1]
                    y_array = [0 -im; im 0]
                    z_array = diagm([1, -1])
                    eigsi = eigvals(kron(i_array, ho_mat))
                    eigsy = eigvals(kron(y_array, ho_mat))
                    eigsz = eigvals(kron(z_array, ho_mat))
                    eigsh = [-70.90875406, -31.04969387, 0, 3.26468993, 38.693758]
                    eigsii = eigvals(kron(i_array, kron(i_array, ho_mat3)))
                    obs_targets = [0, 1, 2]
                    @testset "Obs $obs" for (obs, expected_mean, expected_var, expected_eigs) in
                                            [
                        (Braket.Observables.I() * ho, meani, vari, eigsi),
                        (Braket.Observables.Y() * ho, meany, vary, eigsy),
                        (Braket.Observables.Z() * ho, meanz, varz, eigsz),
                        (ho2 * ho, meanh, varh, eigsh),
                        (
                            Braket.Observables.HermitianObservable(kron(ho_mat2, ho_mat)),
                            meanh,
                            varh,
                            eigsh,
                        ),
                        (Braket.Observables.I() * Braket.Observables.I() * ho3, meanii, varii, eigsii),
                        (Braket.Observables.I() * ho4, meanii, varii, eigsii),
                    ]
                        device = DEVICE
                        shots = SHOTS
                        circuit = three_qubit_circuit(θ, ϕ, φ, obs, obs_targets)
                        shots > 0 && circuit(Braket.Sample, obs, obs_targets)
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        for task in tasks
                            res = result(device(task, shots = shots))
                            variance_expectation_sample_result(
                                res,
                                shots,
                                expected_var,
                                expected_mean,
                                expected_eigs,
                            )
                        end
                    end
                end
                @testset "Result types single Hermitian" begin
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    ho_mat = [
                        -6 2+im -3 -5+2im
                        2-im 0 2-im -5+4im
                        -3 2+im 0 -4+3im
                        -5-2im -5-4im -4-3im -6
                    ]
                    @test ishermitian(ho_mat)
                    ho = Braket.Observables.HermitianObservable(ComplexF64.(ho_mat))
                    ho_mat2 = [1 2; 2 4]
                    ho2 = Braket.Observables.HermitianObservable(ho_mat2)
                    ho_mat3 = [-6 2+im; 2-im 0]
                    ho3 = Braket.Observables.HermitianObservable(ho_mat3)
                    ho_mat4 = kron([1 0; 0 1], [-6 2+im; 2-im 0])
                    ho4 = Braket.Observables.HermitianObservable(ho_mat4)
                    h = Braket.Observables.HermitianObservable(kron(ho_mat2, ho_mat))
                    meani = -5.7267957792059345
                    meanh = -4.30215023196904
                    meanii = -5.78059066879935

                    vari = 43.33800156673375
                    varh = 370.71292282796804
                    varii = 6.268315532585994

                    i_array = [1 0; 0 1]
                    eigsi = eigvals(kron(i_array, ho_mat))
                    eigsh = [-70.90875406, -31.04969387, 0, 3.26468993, 38.693758]
                    eigsii = eigvals(kron(i_array, kron(i_array, ho_mat3)))
                    obs_targets = [0, 1, 2]
                    @testset "Obs $obs" for (
                        obs,
                        targets,
                        expected_mean,
                        expected_var,
                        expected_eigs,
                    ) in [
                        (ho, [1, 2], meani, vari, eigsi),
                        (h, [0, 1, 2], meanh, varh, eigsh),
                        (ho3, [2], meanii, varii, eigsii),
                        (ho4, [1, 2], meanii, varii, eigsii),
                    ]
                        device = DEVICE
                        shots = SHOTS
                        circuit = three_qubit_circuit(θ, ϕ, φ, obs, targets)
                        shots > 0 && circuit(Braket.Sample, obs, targets)
                        tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                        for task in tasks
                            res = result(device(task; shots = shots))
                            variance_expectation_sample_result(
                                res,
                                shots,
                                expected_var,
                                expected_mean,
                                expected_eigs,
                            )
                        end
                    end

                end
                @testset "Result types all selected" begin
                    θ = 0.543
                    ho_mat = [1 2im; -2im 0]
                    ho = Braket.Observables.HermitianObservable(ho_mat)
                    expected_mean = 2 * sin(θ) + 0.5 * cos(θ) + 0.5
                    var_ = 0.25 * (sin(θ) - 4 * cos(θ))^2
                    expected_var = [var_, var_]
                    expected_eigs = eigvals(Hermitian(ho_mat))
                    device = DEVICE
                    shots = SHOTS
                    circuit =
                        Braket.Circuit([(Braket.Rx, 0, θ), (Braket.Rx, 1, θ), (Braket.Variance, ho), (Braket.Expectation, ho, 0)])
                    shots > 0 && circuit(Braket.Sample, ho, 1)
                    tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    for task in tasks 
                        res = result(device(task; shots = shots))
                        tol = get_tol(shots)
                        variance = res.values[1]
                        expectation = res.values[2]
                        if shots > 0
                            samples = res.values[3]
                            @test isapprox(
                                sort(collect(unique(samples))),
                                sort(collect(unique(expected_eigs))),
                                rtol = tol["rtol"],
                                atol = tol["atol"],
                            )
                            @test isapprox(
                                mean(samples),
                                expected_mean,
                                rtol = tol["rtol"],
                                atol = tol["atol"],
                            )
                            @test isapprox(
                                var(samples),
                                var_,
                                rtol = tol["rtol"],
                                atol = tol["atol"],
                            )
                        end
                        @test isapprox(
                            expectation,
                            expected_mean,
                            rtol = tol["rtol"],
                            atol = tol["atol"],
                        )
                        @test isapprox(
                            variance,
                            expected_var,
                            rtol = tol["rtol"],
                            atol = tol["atol"],
                        )
                    end
                end
                @testset "Result types noncommuting" begin
                    shots = 0
                    θ = 0.432
                    ϕ = 0.123
                    φ = -0.543
                    ho_mat = [
                        -6 2+im -3 -5+2im
                        2-im 0 2-im -5+4im
                        -3 2+im 0 -4+3im
                        -5-2im -5-4im -4-3im -6
                    ]
                    obs1 = Braket.Observables.X() * Braket.Observables.Y()
                    obs1_targets = [0, 2]
                    obs2 = Braket.Observables.Z() * Braket.Observables.Z()
                    obs2_targets = [0, 2]
                    obs3 = Braket.Observables.Y() * Braket.Observables.HermitianObservable(ho_mat)
                    obs3_targets = [0, 1, 2]
                    obs3_targets = [0, 1, 2]
                    circuit = three_qubit_circuit(θ, ϕ, φ, obs1, obs1_targets)
                    circuit(Braket.Expectation, obs2, obs2_targets)
                    circuit(Braket.Expectation, obs3, obs3_targets)

                    expected_mean1 = sin(θ) * sin(ϕ) * sin(φ)
                    expected_var1 =
                        (
                            8 * sin(θ)^2 * cos(2φ) * sin(ϕ)^2 - cos(2(θ - ϕ)) - cos(2(θ + ϕ)) +
                            2 * cos(2θ) +
                            2 * cos(2ϕ) +
                            14
                        ) / 16
                    expected_mean2 = 0.849694136476246
                    expected_mean3 = 1.4499810303182408

                    tasks = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    @testset for task in tasks
                        device = DEVICE
                        res = result(device(task, shots = shots))
                        @test isapprox(res.values[1], expected_var1)
                        @test isapprox(res.values[2], expected_mean1)
                        @test isapprox(res.values[3], expected_mean2)
                        @test isapprox(res.values[4], expected_mean3)
                    end
                end
                @testset "Result types noncommuting flipped targets" begin
                    circuit = bell_circ()
                    tp      = Braket.Observables.TensorProduct(["h", "x"])
                    circuit = Braket.Expectation(circuit, tp, [0, 1])
                    circuit = Braket.Expectation(circuit, tp, [1, 0])
                    tasks   = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    @testset for task in tasks
                        device = DEVICE
                        res = result(device(task, shots = 0))
                        @test isapprox(res.values[1], √2 / 2)
                        @test isapprox(res.values[2], √2 / 2)
                    end
                end
                @testset "Result types all noncommuting" begin
                    circuit = bell_circ()
                    ho = [1 2im; -2im 0]
                    circuit(Braket.Expectation, Braket.Observables.HermitianObservable(ho))
                    circuit(Braket.Expectation, Braket.Observables.X())
                    tasks   = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    @testset for task in tasks
                        device = DEVICE
                        res = result(device(task, shots = 0))
                        @test isapprox(res.values[1], [0.5, 0.5])
                        @test isapprox(res.values[2], [0, 0])
                    end
                end
                @testset "Result types observable not in instructions" begin
                    bell = bell_circ()
                    bell(Braket.Expectation, Braket.Observables.X(), 2)
                    bell(Braket.Variance, Braket.Observables.Y(), 3)
                    bell_qasm = ir(bell, Val(:OpenQASM))
                    @test qubit_count(bell) == 4
                    shots = SHOTS
                    device = DEVICE
                    @testset for task in (bell, ir(bell, Val(:JAQCD)), bell_qasm)
                        tol = get_tol(shots)
                        res = result(device(task, shots = shots))
                        @test isapprox(res.values[1], 0, rtol = tol["rtol"], atol = tol["atol"])
                        @test isapprox(res.values[2], 1, rtol = tol["rtol"], atol = tol["atol"])
                    end
                end
            end
        end
        @testset for DEVICE in PURE_DEVICES, SHOTS in SHOT_LIST
            @testset "Result types no shots" begin
                @testset for include_amplitude in [true, false]
                    circuit = bell_circ()
                    circuit(Braket.Expectation, Braket.Observables.H() * Braket.Observables.X(), 0, 1)
                    include_amplitude && circuit(Braket.Amplitude, ["01", "10", "00", "11"])
                    tasks   = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    @testset for task in tasks
                        device = DEVICE
                        shots = 0
                        res = result(device(task, shots = 0))
                        @test length(res.result_types) == (include_amplitude ? 2 : 1)
                        @test isapprox(
                            res[Braket.Expectation(Braket.Observables.H() * Braket.Observables.X(), [0, 1])],
                            1 / √2,
                        )
                        if include_amplitude
                            amps = res[Braket.Amplitude(["01", "10", "00", "11"])]
                            @test isapprox(amps["01"], 0)
                            @test isapprox(amps["10"], 0)
                            @test isapprox(amps["00"], 1 / √2)
                            @test isapprox(amps["11"], 1 / √2)
                        end
                    end
                end
            end
            if SHOTS > 0
                @testset "Multithreaded Bell pair" begin
                    tol = get_tol(SHOTS)
                    tasks = (bell_circ, (()->ir(bell_circ(), Val(:JAQCD))), (()->ir(bell_circ(), Val(:OpenQASM))))
                    device = DEVICE
                    @testset for task in tasks, task_count in (1, 10)
                        task_array = [task() for ii = 1:task_count]
                        batch_results = results(device(task_array, shots=SHOTS))
                        for r in batch_results
                            @test isapprox(
                                r.measurement_probabilities["00"],
                                0.5,
                                rtol = tol["rtol"],
                                atol = tol["atol"],
                            )
                            @test isapprox(
                                r.measurement_probabilities["11"],
                                0.5,
                                rtol = tol["rtol"],
                                atol = tol["atol"],
                            )
                            @test length(r.measurements) == SHOTS
                        end
                    end
                end
            end
        end
        @testset for DEVICE in NOISE_DEVICES, SHOTS in SHOT_LIST
            @testset "noisy circuit 1 qubit noise full probability" begin
                shots = SHOTS
                tol = get_tol(shots)
                circuit = Braket.Circuit([(Braket.X, 0), (Braket.X, 1), (Braket.BitFlip, 0, 0.1), (Braket.Probability,)])
                tasks   = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                device  = DEVICE
                for task in tasks
                    res = result(device(task, shots = shots))
                    @test length(res.result_types) == 1
                    @test isapprox(
                        res[Probability()],
                        [0.0, 0.1, 0, 0.9],
                        rtol = tol["rtol"],
                        atol = tol["atol"],
                    )
                end
            end
            @testset "noisy circuit 2 qubit noise full probability" begin
                shots = SHOTS
                tol = get_tol(shots)
                K0 = √0.9 * diagm(ones(4))
                K1 = √0.1 * kron([0.0 1.0; 1.0 0.0], [0.0 1.0; 1.0 0.0])
                circuit =
                    Braket.Circuit([(Braket.X, 0), (Braket.X, 1), (Braket.Kraus, [0, 1], [K0, K1]), (Braket.Probability,)])
                tasks   = (circuit, ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                device = DEVICE
                for task in tasks
                    res = result(device(task, shots = shots))
                    @test length(res.result_types) == 1
                    @test isapprox(
                        res[Braket.Probability()],
                        [0.1, 0.0, 0, 0.9],
                        rtol = tol["rtol"],
                        atol = tol["atol"],
                    )
                end
            end
        end
    end
end
