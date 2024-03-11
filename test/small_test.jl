using Test,
    cuStateVec,
    CUDA,
    Statistics,
    LinearAlgebra,
    Braket,
    Braket.Observables,
    BraketStateVector
using Braket: I, name

@testset "Correctness" begin
    Braket.IRType[] = :JAQCD
    PURE_DEVICE = LocalSimulator("braket_sv")
    NOISE_DEVICE = LocalSimulator("braket_dm")
    ALL_DEVICES = [PURE_DEVICE, NOISE_DEVICE]
    PURE_DEVICES = [PURE_DEVICE]
    NOISE_DEVICES = [NOISE_DEVICE]
    if CUDA.functional()
        CU_PURE_DEVICE = LocalSimulator("braket_sv_custatevec")
        CU_NOISE_DEVICE = LocalSimulator("braket_dm_custatevec")
        append!(ALL_DEVICES, [CU_PURE_DEVICE, CU_NOISE_DEVICE])
        push!(PURE_DEVICES, CU_PURE_DEVICE)
        push!(NOISE_DEVICES, CU_NOISE_DEVICE)
    end
    SHOT_LIST = (8000,)

    get_tol(shots::Int) = return (
        shots > 0 ? Dict("atol" => 0.1, "rtol" => 0.15) : Dict("atol" => 0.01, "rtol" => 0)
    )

    bell_circ() = Circuit([(H, 0), (CNot, 0, 1)])
    three_qubit_circuit(
        θ::Float64,
        ϕ::Float64,
        φ::Float64,
        obs::Observables.Observable,
        obs_targets::Vector{Int},
    ) = Circuit([
        (Rx, 0, θ),
        (Rx, 1, ϕ),
        (Rx, 2, φ),
        (CNot, 0, 1),
        (CNot, 1, 2),
        (Variance, obs, obs_targets),
        (Expectation, obs, obs_targets),
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
            ("default", "StateVectorSimulator"),
            ("braket_sv", "StateVectorSimulator"),
            ("braket_dm", "DensityMatrixSimulator"),
        ]
            local_simulator_device = LocalSimulator(backend)
            @test name(local_simulator_device) == device_name
        end
        @testset "Device $DEVICE, shots $SHOTS" for DEVICE in ALL_DEVICES,
            SHOTS in SHOT_LIST
            @testset "Result types all selected" begin
                θ = 0.543
                ho_mat = [1 2im; -2im 0]
                ho = Observables.HermitianObservable(ho_mat)
                expected_mean = 2 * sin(θ) + 0.5 * cos(θ) + 0.5
                var_ = 0.25 * (sin(θ) - 4 * cos(θ))^2
                expected_var = [var_, var_]
                expected_eigs = eigvals(Hermitian(ho_mat))
                device = DEVICE
                shots = SHOTS
                circuit =
                    Circuit([(Rx, 0, θ), (Rx, 1, θ), (Variance, ho), (Expectation, ho, 0)])
                #shots > 0 && circuit(Sample, ho, 1)
                @testset for task in (ir(circuit, Val(:JAQCD)), ir(circuit, Val(:OpenQASM)))
                    @show task
                    res         = result(device(task; shots = shots))
                    tol         = get_tol(shots)
                    variance    = res.values[1]
                    expectation = res.values[2]
                    #=if shots > 0
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
                    end=#
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
        end
    end
end
