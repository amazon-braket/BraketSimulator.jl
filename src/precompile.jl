function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{HermitianObservable, Vector{Int64}}}})   # time: 2.3760564
    Base.precompile(Tuple{var"#556#threadsfor_fun#215"{var"#556#threadsfor_fun#211#216"{Matrix{ComplexF64}, Int64, Vector{Vector{Int64}}, Vector{Int64}, Transpose{ComplexF64, Matrix{ComplexF64}}, UnitRange{Int64}}},Int64})   # time: 0.14115863
    Base.precompile(Tuple{var"#357#threadsfor_fun#123"{var"#357#threadsfor_fun#120#124"{16, ComplexF64, SMatrix{16, 16, ComplexF64, 256}, Vector{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.12406683
    Base.precompile(Tuple{var"#403#threadsfor_fun#143"{var"#403#threadsfor_fun#140#144"{Matrix{ComplexF64}, Int64, Int64, Int64, Int64, Tuple{SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}}, Tuple{SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}}, UnitRange{Int64}}},Int64})   # time: 0.11665836
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{HermitianObservable, Vector{Int64}}}})   # time: 0.08039276
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},DoubleExcitation,Int64,Vararg{Int64}})   # time: 0.0797019
    Base.precompile(Tuple{var"#536#threadsfor_fun#204"{var"#536#threadsfor_fun#201#205"{Vector{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, Matrix{ComplexF64}, UnitRange{Int64}}},Int64})   # time: 0.07728745
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Kraus,Int64,Vararg{Int64}})   # time: 0.06592051
    Base.precompile(Tuple{var"#403#threadsfor_fun#143"{var"#403#threadsfor_fun#140#144"{Matrix{ComplexF64}, Int64, Int64, Int64, Int64, Tuple{SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}}, Tuple{SMatrix{4, 4, ComplexF64, 16}, SMatrix{4, 4, ComplexF64, 16}}, UnitRange{Int64}}},Int64})   # time: 0.062411573
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots,), Tuple{Int64}},typeof(simulate),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},OpenQasmProgram})   # time: 0.06121916
    Base.precompile(Tuple{var"#403#threadsfor_fun#143"{var"#403#threadsfor_fun#140#144"{Matrix{ComplexF64}, Int64, Int64, Int64, Int64, NTuple{16, SMatrix{4, 4, ComplexF64, 16}}, NTuple{16, SMatrix{4, 4, ComplexF64, 16}}, UnitRange{Int64}}},Int64})   # time: 0.052433595
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.X, Int64}}})   # time: 0.050354417
    Base.precompile(Tuple{typeof(calculate),Braket.DensityMatrix,StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.047030203
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Int64})   # time: 0.045860793
    Base.precompile(Tuple{typeof(calculate),Amplitude,StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.04560804
    isdefined(BraketSimulator, Symbol("#60#64")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#60#64")),Instruction{BitFlip}})   # time: 0.045065258
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{HermitianObservable},Int64,Int64,Vararg{Int64}})   # time: 0.044758502
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TwoQubitDepolarizing,Int64,Vararg{Int64}})   # time: 0.044507124
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.StandardObservable}, Tuple{Int64, Int64}}}})   # time: 0.04418025
    Base.precompile(Tuple{typeof(calculate),Probability,StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.041471716
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.NonCompositeObservable}, Vector{Int64}}}})   # time: 0.0414185
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.StandardObservable}, Vector{Int64}}}})   # time: 0.04043337
    Base.precompile(Tuple{typeof(^),CCNot,Int64})   # time: 0.03875475
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiQubitPhaseShift{2},Int64,Vararg{Int64}})   # time: 0.03686717
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.StandardObservable, Vector{Int64}}}})   # time: 0.03636653
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiQubitPauliChannel{2},Int64,Vararg{Int64}})   # time: 0.035675667
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots,), Tuple{Int64}},typeof(simulate),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},OpenQasmProgram})   # time: 0.03539539
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.NonCompositeObservable}, NTuple{4, Int64}}}})   # time: 0.034470204
    isdefined(BraketSimulator, Symbol("#162#163")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#162#163")),Tuple{Braket.Observables.Y, Int64}})   # time: 0.03252417
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),D<:BraketSimulator.AbstractSimulator,Union{Circuit, Program},Val{:OpenQASM}})   # time: 0.031343497
    Base.precompile(Tuple{typeof(apply_gate!),DoubleExcitation,Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.030710414
    Base.precompile(Tuple{typeof(apply_gate!),Unitary,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.030677125
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiQubitPhaseShift{4},Int64,Vararg{Int64}})   # time: 0.029619085
    Base.precompile(Tuple{var"#403#threadsfor_fun#143"{var"#403#threadsfor_fun#140#144"{Matrix{ComplexF64}, Int64, Int64, Int64, Int64, NTuple{4, SMatrix{4, 4, ComplexF64, 16}}, NTuple{4, SMatrix{4, 4, ComplexF64, 16}}, UnitRange{Int64}}},Int64})   # time: 0.028585028
    Base.precompile(Tuple{typeof(matrix_rep),DoubleExcitation})   # time: 0.027744127
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},SingleExcitation,Int64,Vararg{Int64}})   # time: 0.026039707
    Base.precompile(Tuple{var"#388#threadsfor_fun#136"{var"#388#threadsfor_fun#133#137"{Matrix{ComplexF64}, Int64, NTuple{4, SMatrix{2, 2, ComplexF64, 4}}, NTuple{4, SMatrix{2, 2, ComplexF64, 4}}, UnitRange{Int64}}},Int64})   # time: 0.025726128
    Base.precompile(Tuple{typeof(apply_gate!),MultiQubitPhaseShift{4},Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.025318487
    Base.precompile(Tuple{typeof(apply_gate!),MultiRZ,Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.024543203
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.StandardObservable}, Tuple{Int64, Int64}}}})   # time: 0.024210926
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.NonCompositeObservable}, Tuple{Int64, Int64}}}})   # time: 0.02415638
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{HermitianObservable},Int64,Int64,Vararg{Int64}})   # time: 0.023168875
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Unitary,Int64,Vararg{Int64}})   # time: 0.022914968
    Base.precompile(Tuple{typeof(apply_gate!),Control{X, 2},Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.022110343
    Base.precompile(Tuple{var"#357#threadsfor_fun#123"{var"#357#threadsfor_fun#120#124"{8, ComplexF64, SMatrix{8, 8, ComplexF64, 64}, Vector{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.022048157
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Int64}}})   # time: 0.021454116
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Braket.DensityMatrix},DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.021367287
    Base.precompile(Tuple{typeof(apply_gate!),MultiQubitPhaseShift{3},Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.0211587
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.H, Int64}}})   # time: 0.021102335
    Base.precompile(Tuple{var"#357#threadsfor_fun#123"{var"#357#threadsfor_fun#120#124"{8, ComplexF64, Diagonal{ComplexF64, SVector{8, ComplexF64}}, Vector{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.02088434
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiQubitPhaseShift{3},Int64,Vararg{Int64}})   # time: 0.020738542
    Base.precompile(Tuple{typeof(apply_gate!),Val{true},MultiRZ,Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.018903581
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),D<:BraketSimulator.AbstractSimulator,Union{Circuit, Program},Val{:JAQCD}})   # time: 0.01855247
    Base.precompile(Tuple{typeof(samples),StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.018172288
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.X, Vector{Int64}}}})   # time: 0.018068293
    Base.precompile(Tuple{typeof(apply_gate!),MultiRZ,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.01697871
    Base.precompile(Tuple{var"#357#threadsfor_fun#123"{var"#357#threadsfor_fun#120#124"{16, ComplexF64, Diagonal{ComplexF64, SVector{16, ComplexF64}}, Vector{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.016731136
    Base.precompile(Tuple{typeof(apply_gate!),MultiQubitPhaseShift{2},Vector{ComplexF64},Int64,Int64})   # time: 0.015248993
    Base.precompile(Tuple{typeof(apply_gate!),Val{true},Unitary,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.014217445
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Tuple{Int64}}}})   # time: 0.0140325
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.X, Vector{Int64}}}})   # time: 0.014004589
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.Z, Vector{Int64}}}})   # time: 0.013751209
    Base.precompile(Tuple{typeof(_validate_ir_results_compatibility),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{NamedTuple{(:type,), Tuple{String}}},Val{:OpenQASM}})   # time: 0.013238872
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},CNot,Int64,Vararg{Int64}})   # time: 0.013133743
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Tuple{Int64}}}})   # time: 0.012980301
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.H, Vector{Int64}}}})   # time: 0.012828214
    Base.precompile(Tuple{var"#419#threadsfor_fun#151"{var"#419#threadsfor_fun#147#152"{Matrix{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, Tuple{Adjoint{ComplexF64, Matrix{ComplexF64}}}, Vector{Matrix{ComplexF64}}, UnitRange{Int64}}},Int64})   # time: 0.012743572
    Base.precompile(Tuple{var"#506#threadsfor_fun#195"{var"#506#threadsfor_fun#194#196"{Vector{ComplexF64}, Int64, Matrix{ComplexF64}, UnitRange{Int64}}},Int64})   # time: 0.012427671
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.H, Vector{Int64}}}})   # time: 0.012073245
    Base.precompile(Tuple{typeof(_validate_ir_results_compatibility),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{NamedTuple{(:type,), Tuple{String}}},Val{:JAQCD}})   # time: 0.011608869
    Base.precompile(Tuple{typeof(apply_gate!),Unitary,Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.011561283
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Int64}}})   # time: 0.011552837
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.H, Int64}}})   # time: 0.011530246
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.X, Int64}}})   # time: 0.01139113
    Base.precompile(Tuple{typeof(expectation),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},HermitianObservable,Int64})   # time: 0.011389514
    Base.precompile(Tuple{typeof(matrix_rep),MultiQubitPhaseShift{4}})   # time: 0.011219284
    Base.precompile(Tuple{typeof(apply_gate!),Val{true},MultiRZ,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.010829912
    Base.precompile(Tuple{var"#441#threadsfor_fun#171"{var"#441#threadsfor_fun#170#172"{Vector{Float64}, Vector{Vector{Int64}}, Vector{Float64}, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.010822869
    Base.precompile(Tuple{var"#220#threadsfor_fun#75"{var"#220#threadsfor_fun#74#76"{ComplexF64, SMatrix{4, 4, ComplexF64, 16}, Vector{ComplexF64}, Vector{UnitRange{Int64}}, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.009974208
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.Z, Vector{Int64}}}})   # time: 0.009941416
    Base.precompile(Tuple{var"#521#threadsfor_fun#198"{var"#521#threadsfor_fun#197#199"{Vector{ComplexF64}, Int64, Int64, Matrix{ComplexF64}, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.009609202
    Base.precompile(Tuple{typeof(apply_noise!),Kraus,Matrix{ComplexF64},Int64,Int64,Int64,Int64,Int64})   # time: 0.009606712
    Base.precompile(Tuple{var"#321#threadsfor_fun#91"{var"#321#threadsfor_fun#87#92"{ComplexF64, Vector{ComplexF64}, Bool, SMatrix{4, 4, ComplexF64, 16}, Int64, Int64, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.009529764
    Base.precompile(Tuple{var"#472#threadsfor_fun#181"{var"#472#threadsfor_fun#178#182"{Matrix{ComplexF64}, Int64, Vector{Int64}, Matrix{ComplexF64}, UnitRange{Int64}}},Int64})   # time: 0.009195421
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots, :inputs), Tuple{Int64, Dict{String, Float64}}},typeof(simulate),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Program,Int64})   # time: 0.009151381
    Base.precompile(Tuple{typeof(apply_gate!),CNot,Vector{ComplexF64},Int64,Int64})   # time: 0.009115913
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots, :inputs), Tuple{Int64, Dict{String, Float64}}},typeof(simulate),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Program,Int64})   # time: 0.008724961
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Result},StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.008643294
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots, :inputs), Tuple{Int64, Dict{String, Float64}}},typeof(simulate),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},OpenQasmProgram})   # time: 0.008624748
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Vector{Int64}}}})   # time: 0.008502583
    Base.precompile(Tuple{typeof(calculate),Expectation,StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.008263838
    Base.precompile(Tuple{typeof(^),CNot,Int64})   # time: 0.008205919
    Base.precompile(Tuple{typeof(Core.kwcall),NamedTuple{(:shots, :inputs), Tuple{Int64, Dict{String, Float64}}},typeof(simulate),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},OpenQasmProgram})   # time: 0.008148667
    Base.precompile(Tuple{typeof(apply_gate!),H,Vector{ComplexF64},Int64})   # time: 0.008035914
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},BitFlip,Int64})   # time: 0.007911752
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Int64})   # time: 0.00772037
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},HermitianObservable,Int64,Int64,Vararg{Int64}})   # time: 0.007329834
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{HermitianObservable, Int64}}})   # time: 0.00700907
    Base.precompile(Tuple{var"#457#threadsfor_fun#175"{var"#457#threadsfor_fun#173#176"{Vector{Float64}, Int64, Vector{Int64}, Vector{Float64}, UnitRange{Int64}}},Int64})   # time: 0.00699316
    Base.precompile(Tuple{typeof(expectation),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},HermitianObservable,Int64,Vararg{Int64}})   # time: 0.006969504
    Base.precompile(Tuple{typeof(matrix_rep),MultiQubitPhaseShift{3}})   # time: 0.006872365
    Base.precompile(Tuple{typeof(marginal_probability),Vector{Float64},Int64,Int64})   # time: 0.006763912
    Base.precompile(Tuple{typeof(marginal_probability),Vector{Float64},Int64,Tuple{Int64, Int64}})   # time: 0.006747198
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.I, Vector{Int64}}}})   # time: 0.006616665
    Base.precompile(Tuple{Type{StateVectorSimulator},Int64,Int64})   # time: 0.00658437
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.Z, Int64}}})   # time: 0.006555793
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.I, Int64}}})   # time: 0.006480999
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CSwap}}})   # time: 0.006389134
    Base.precompile(Tuple{typeof(marginal_probability),Vector{Float64},Int64,Tuple{Int64}})   # time: 0.006313712
    Base.precompile(Tuple{typeof(probabilities),StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.006088001
    Base.precompile(Tuple{typeof(apply_observables!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Tuple{Braket.Observables.StandardObservable, Vector{Int64}}}})   # time: 0.006057755
    Base.precompile(Tuple{var"#220#threadsfor_fun#75"{var"#220#threadsfor_fun#74#76"{ComplexF64, Diagonal{ComplexF64, SVector{4, ComplexF64}}, Vector{ComplexF64}, Vector{UnitRange{Int64}}, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.005851834
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},HermitianObservable,Int64,Int64})   # time: 0.005848543
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiRZ,Int64,Vararg{Int64}})   # time: 0.005818084
    Base.precompile(Tuple{typeof(apply_gate!),MultiQubitPhaseShift{1},Vector{ComplexF64},Int64})   # time: 0.005786125
    isdefined(BraketSimulator, Symbol("#60#64")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#60#64")),Instruction{Rx}})   # time: 0.00568604
    Base.precompile(Tuple{var"#388#threadsfor_fun#136"{var"#388#threadsfor_fun#133#137"{Matrix{ComplexF64}, Int64, Tuple{SMatrix{2, 2, ComplexF64, 4}, SMatrix{2, 2, ComplexF64, 4}}, Tuple{SMatrix{2, 2, ComplexF64, 4}, SMatrix{2, 2, ComplexF64, 4}}, UnitRange{Int64}}},Int64})   # time: 0.005640168
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{HermitianObservable, Int64}}})   # time: 0.005586166
    Base.precompile(Tuple{typeof(probabilities),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.005552664
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.Z},Int64,Int64})   # time: 0.005303585
    Base.precompile(Tuple{typeof(evolve!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Instruction}})   # time: 0.005119041
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CZ}}})   # time: 0.005095164
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},H,Int64})   # time: 0.005078246
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CY}}})   # time: 0.005051496
    isdefined(BraketSimulator, Symbol("#60#64")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#60#64")),Instruction{Ry}})   # time: 0.004877292
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Control{X, 1},Int64,Vararg{Int64}})   # time: 0.004856696
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},CSwap,Int64,Vararg{Int64}})   # time: 0.004838249
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.I, Vector{Int64}}}})   # time: 0.004632425
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.Z, Int64}}})   # time: 0.004622792
    Base.precompile(Tuple{typeof(matrix_rep),MultiQubitPhaseShift{2}})   # time: 0.00461696
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.I, Int64}}})   # time: 0.004577912
    Base.precompile(Tuple{var"#336#threadsfor_fun#93"{var"#336#threadsfor_fun#88#94"{Vector{ComplexF64}, Tuple{Int64, Int64}, ComplexF64, ComplexF64, ComplexF64, ComplexF64, Int64, Int64, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.004202628
    Base.precompile(Tuple{typeof(matrix_rep),Control{MultiQubitPhaseShift{2}, 2}})   # time: 0.004102118
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{Braket.Observables.Y, Vector{Int64}}}})   # time: 0.004091999
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},MultiQubitPhaseShift{1},Int64})   # time: 0.004060454
    Base.precompile(Tuple{typeof(evolve!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Instruction{Kraus}}})   # time: 0.004041914
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Circuit,Val{:OpenQASM}})   # time: 0.004018167
    Base.precompile(Tuple{typeof(matrix_rep),SingleExcitation})   # time: 0.004016366
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Probability},DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.003991788
    Base.precompile(Tuple{typeof(expectation),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{HermitianObservable},Int64,Vararg{Int64}})   # time: 0.00398108
    Base.precompile(Tuple{typeof(calculate),Expectation,DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.003972124
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Int64,Vararg{Int64}})   # time: 0.003953711
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Result},DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.00392808
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},CY,Int64,Vararg{Int64}})   # time: 0.003895222
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},CZ,Int64,Vararg{Int64}})   # time: 0.003787091
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},HermitianObservable,Int64,Int64})   # time: 0.003721872
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Probability},StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.003683374
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Circuit,Val{:JAQCD}})   # time: 0.00366429
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Int64,Vararg{Int64}})   # time: 0.003664
    Base.precompile(Tuple{var"#251#threadsfor_fun#82"{var"#251#threadsfor_fun#78#83"{ComplexF64, Vector{ComplexF64}, Bool, SMatrix{4, 4, ComplexF64, 16}, Int64, Int64, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.003619081
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Expectation},DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.003612213
    Base.precompile(Tuple{typeof(apply_observables!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Tuple{TensorProduct{Braket.Observables.NonCompositeObservable}, Tuple{Int64, Int64}}}})   # time: 0.003567337
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{HermitianObservable},Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.003555792
    Base.precompile(Tuple{typeof(_generate_results),Vector{Braket.IR.AbstractProgramResult},Vector{Expectation},StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.003515872
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Control{X, 2},Int64,Vararg{Int64}})   # time: 0.003494703
    Base.precompile(Tuple{typeof(apply_gate!),MultiRZ,Vector{ComplexF64},Int64,Int64})   # time: 0.003429628
    Base.precompile(Tuple{typeof(matrix_rep),MultiQubitPhaseShift{1}})   # time: 0.003388248
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Int64})   # time: 0.003384418
    Base.precompile(Tuple{typeof(calculate),Variance,DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.003176043
    Base.precompile(Tuple{typeof(calculate),Variance,StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.003117915
    isdefined(BraketSimulator, Symbol("#162#163")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#162#163")),Tuple{Braket.Observables.H, Int64}})   # time: 0.003023545
    Base.precompile(Tuple{typeof(_validate_input_provided),Circuit})   # time: 0.003010547
    isdefined(BraketSimulator, Symbol("#207#208")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#207#208")),Tuple{Braket.Observables.H, Vector{Int64}}})   # time: 0.003008123
    isdefined(BraketSimulator, Symbol("#162#163")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#162#163")),Tuple{Braket.Observables.X, Int64}})   # time: 0.002941127
    Base.precompile(Tuple{typeof(apply_gate!),Rz,Vector{ComplexF64},Int64})   # time: 0.002940584
    Base.precompile(Tuple{var"#419#threadsfor_fun#151"{var"#419#threadsfor_fun#147#152"{Matrix{ComplexF64}, Vector{Vector{Int64}}, Vector{Int64}, Tuple{Adjoint{ComplexF64, Matrix{ComplexF64}}, Adjoint{ComplexF64, Matrix{ComplexF64}}}, Vector{Matrix{ComplexF64}}, UnitRange{Int64}}},Int64})   # time: 0.002931998
    Base.precompile(Tuple{typeof(calculate),Variance,AbstractSimulator})   # time: 0.002875833
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},CCNot,Int64,Vararg{Int64}})   # time: 0.002809754
    Base.precompile(Tuple{typeof(^),PSwap,Int64})   # time: 0.002741585
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Unitary}}})   # time: 0.002696421
    Base.precompile(Tuple{typeof(apply_gate!),CCNot,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.002666251
    Base.precompile(Tuple{typeof(apply_observable!),HermitianObservable,Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.002656832
    Base.precompile(Tuple{typeof(expectation),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Vararg{Int64}})   # time: 0.002642792
    Base.precompile(Tuple{typeof(expectation),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Vararg{Int64}})   # time: 0.002624209
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction}})   # time: 0.002521376
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Control{MultiQubitPhaseShift{2}, 2},Int64,Vararg{Int64}})   # time: 0.002457497
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{YY}}})   # time: 0.002411541
    Base.precompile(Tuple{typeof(matrix_rep),Control{MultiQubitPhaseShift{1}, 1}})   # time: 0.002357129
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{HermitianObservable},Matrix{ComplexF64},Int64,Int64,Int64})   # time: 0.002340128
    Base.precompile(Tuple{typeof(expectation),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.NonCompositeObservable},Int64,Vararg{Int64}})   # time: 0.002331834
    Base.precompile(Tuple{typeof(expectation),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Vararg{Int64}})   # time: 0.002299791
    Base.precompile(Tuple{typeof(apply_observable!),HermitianObservable,Matrix{ComplexF64},Int64,Int64,Int64})   # time: 0.0022985
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},HermitianObservable,Int64,Int64,Vararg{Int64}})   # time: 0.00229296
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{H}}})   # time: 0.002282665
    Base.precompile(Tuple{typeof(apply_gate!),Control{X, 1},Vector{ComplexF64},Int64,Int64})   # time: 0.002236294
    Base.precompile(Tuple{typeof(apply_gate!),Control{I, 1},Vector{ComplexF64},Int64,Int64})   # time: 0.00222787
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},HermitianObservable,Int64})   # time: 0.002189707
    Base.precompile(Tuple{typeof(^),X,Int64})   # time: 0.002178541
    Base.precompile(Tuple{typeof(apply_gate!),Val{true},Unitary,Vector{ComplexF64},Int64,Int64,Int64,Int64})   # time: 0.002167166
    Base.precompile(Tuple{typeof(matrix_rep),XY})   # time: 0.002118794
    Base.precompile(Tuple{typeof(matrix_rep),XX})   # time: 0.002106209
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{ECR}}})   # time: 0.002103711
    isdefined(BraketSimulator, Symbol("#188#190")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#188#190")),HermitianObservable})   # time: 0.002097624
    Base.precompile(Tuple{typeof(apply_gate!),CV,Vector{ComplexF64},Int64,Int64})   # time: 0.00209163
    Base.precompile(Tuple{typeof(matrix_rep),YY})   # time: 0.002058165
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CCNot}}})   # time: 0.002057918
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{V}}})   # time: 0.002049661
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Vi}}})   # time: 0.002040998
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},PhaseFlip,Int64})   # time: 0.002026789
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{ZZ}}})   # time: 0.002024126
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{XY}}})   # time: 0.002018045
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Int64,Vararg{Int64}})   # time: 0.002012294
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{PSwap}}})   # time: 0.002007076
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},HermitianObservable,Int64})   # time: 0.002006084
    Base.precompile(Tuple{var"#266#threadsfor_fun#84"{var"#266#threadsfor_fun#79#85"{Vector{ComplexF64}, Tuple{Int64, Int64}, ComplexF64, ComplexF64, ComplexF64, ComplexF64, Int64, Int64, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.001982498
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CPhaseShift10}}})   # time: 0.001980336
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{MS}}})   # time: 0.001973913
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CPhaseShift}}})   # time: 0.001968584
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CPhaseShift00}}})   # time: 0.001967121
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Z}}})   # time: 0.001953294
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{I}}})   # time: 0.001952793
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CNot}}})   # time: 0.001940213
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{U}}})   # time: 0.00193508
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Y}}})   # time: 0.001934711
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{CPhaseShift01}}})   # time: 0.001930041
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Swap}}})   # time: 0.001927331
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{GPi}}})   # time: 0.001898796
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{XX}}})   # time: 0.001895626
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{GPi2}}})   # time: 0.001883828
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Rz}}})   # time: 0.001863544
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{X}}})   # time: 0.001833085
    Base.precompile(Tuple{var"#53#threadsfor_fun#40"{var"#53#threadsfor_fun#39#41"{Vector{ComplexF64}, Int64, Int64, Int64, Int64, Float64, Float64, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.001826791
    Base.precompile(Tuple{typeof(permute_probability),Vector{Float64},Int64,Vector{Int64}})   # time: 0.001820791
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{PhaseShift}}})   # time: 0.001816795
    Base.precompile(Tuple{typeof(expectation_op_squared),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Int64})   # time: 0.001796248
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Ry}}})   # time: 0.001789662
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Circuit,Val{:OpenQASM}})   # time: 0.001764041
    Base.precompile(Tuple{var"#204#threadsfor_fun#72"{var"#204#threadsfor_fun#71#73"{Vector{ComplexF64}, ComplexF64, ComplexF64, ComplexF64, ComplexF64, Bool, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.001741707
    Base.precompile(Tuple{typeof(apply_gate!),Control{MultiQubitPhaseShift{2}, 2},Vector{ComplexF64},Int64,Int64})   # time: 0.001721957
    Base.precompile(Tuple{typeof(expectation_op_squared),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},TensorProduct{Braket.Observables.StandardObservable},Int64,Int64,Vararg{Int64}})   # time: 0.001721418
    isdefined(BraketSimulator, Symbol("#51#54")) && Base.precompile(Tuple{getfield(BraketSimulator, Symbol("#51#54")),Tuple{Braket.IR.AbstractProgramResult, Any}})   # time: 0.001702045
    Base.precompile(Tuple{typeof(evolve!),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Vector{Instruction{Rx}}})   # time: 0.001694167
    Base.precompile(Tuple{typeof(samples),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.00169363
    Base.precompile(Tuple{typeof(matrix_rep),Rx})   # time: 0.001682583
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},TwoQubitDephasing,Int64,Vararg{Int64}})   # time: 0.001661709
    Base.precompile(Tuple{typeof(inv),MS})   # time: 0.001660417
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Depolarizing,Int64})   # time: 0.00160096
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.NonCompositeObservable},Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.001573957
    Base.precompile(Tuple{typeof(index_to_endian_bits),Int64,Int64})   # time: 0.001559832
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.StandardObservable},Matrix{ComplexF64},Int64,Int64,Int64})   # time: 0.0015475
    Base.precompile(Tuple{typeof(apply_gate!),Ry,Vector{ComplexF64},Int64})   # time: 0.001546374
    Base.precompile(Tuple{typeof(matrix_rep),I})   # time: 0.001544544
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Control{MultiQubitPhaseShift{1}, 1},Int64})   # time: 0.001539502
    Base.precompile(Tuple{typeof(_validate_operation_qubits),Vector{Instruction{Rx}}})   # time: 0.00153604
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.NonCompositeObservable},Matrix{ComplexF64},Int64,Int64,Int64})   # time: 0.001534041
    Base.precompile(Tuple{var"#306#threadsfor_fun#89"{var"#306#threadsfor_fun#86#90"{Vector{ComplexF64}, ComplexF64, ComplexF64, ComplexF64, ComplexF64, Bool, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.001509835
    Base.precompile(Tuple{typeof(apply_gate!),Control{MultiQubitPhaseShift{1}, 1},Vector{ComplexF64},Int64})   # time: 0.001500166
    Base.precompile(Tuple{typeof(apply_gate!),Val{true},Control{MultiQubitPhaseShift{2}, 2},Vector{ComplexF64},Int64,Int64})   # time: 0.00148471
    Base.precompile(Tuple{typeof(matrix_rep),PhaseShift})   # time: 0.001484083
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.StandardObservable},Vector{ComplexF64},Int64,Int64,Int64})   # time: 0.001468376
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.NonCompositeObservable},Matrix{ComplexF64},Int64,Int64})   # time: 0.001460708
    Base.precompile(Tuple{typeof(expectation),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},HermitianObservable,Int64})   # time: 0.001454329
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.StandardObservable},Matrix{ComplexF64},Int64,Int64})   # time: 0.001452959
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},GeneralizedAmplitudeDamping,Int64})   # time: 0.001450046
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.NonCompositeObservable},Vector{ComplexF64},Int64,Int64})   # time: 0.001447251
    Base.precompile(Tuple{typeof(_validate_ir_instructions_compatibility),StateVectorSimulator{ComplexF64, Vector{ComplexF64}},Circuit,Val{:JAQCD}})   # time: 0.001447165
    Base.precompile(Tuple{typeof(apply_observable!),TensorProduct{Braket.Observables.StandardObservable},Vector{ComplexF64},Int64,Int64})   # time: 0.001421333
    Base.precompile(Tuple{typeof(^),MS,Int64})   # time: 0.001409875
    Base.precompile(Tuple{typeof(matrix_rep),Control{I, 1}})   # time: 0.001371751
    Base.precompile(Tuple{typeof(matrix_rep),Control{X, 1}})   # time: 0.001291998
    Base.precompile(Tuple{var"#373#threadsfor_fun#127"{var"#373#threadsfor_fun#126#128"{PhaseFlip, Matrix{ComplexF64}, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.00128375
    Base.precompile(Tuple{typeof(copy),StateVectorSimulator{ComplexF64, Vector{ComplexF64}}})   # time: 0.001264997
    Base.precompile(Tuple{typeof(calculate),Probability,DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}}})   # time: 0.001251415
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},AmplitudeDamping,Int64})   # time: 0.001234916
    Base.precompile(Tuple{var"#73#threadsfor_fun#43"{var"#73#threadsfor_fun#42#44"{Vector{ComplexF64}, Int64, Int64, Int64, Int64, Float64, Float64, Vector{Int64}, UnitRange{Int64}}},Int64})   # time: 0.001232541
    Base.precompile(Tuple{var"#236#threadsfor_fun#80"{var"#236#threadsfor_fun#77#81"{Vector{ComplexF64}, ComplexF64, ComplexF64, ComplexF64, ComplexF64, Bool, Int64, Int64, Int64, Int64, UnitRange{Int64}}},Int64})   # time: 0.001212125
    Base.precompile(Tuple{typeof(_evolve_op!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},PauliChannel,Int64})   # time: 0.001116707
    Base.precompile(Tuple{Type{DensityMatrixSimulator},Int64,Int64})   # time: 0.001066707
    Base.precompile(Tuple{typeof(evolve!),DensityMatrixSimulator{ComplexF64, Matrix{ComplexF64}},Vector{Instruction{XX}}})   # time: 0.001022039
end
