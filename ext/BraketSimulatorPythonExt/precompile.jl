function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(jl_convert),Type{OpenQasmProgram},Py})   # time: 0.5335359
    Base.precompile(Tuple{StateVectorSimulator{ComplexF64, Vector{ComplexF64}},PyList{Any},Int64})   # time: 0.014752082
    Base.precompile(Tuple{typeof(union_convert),Type,Py})   # time: 0.00300354
end
