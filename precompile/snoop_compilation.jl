using SnoopCompile

inf_timing = @snoopi_deep include(joinpath(@__DIR__, "..", "test", "runtests.jl"))
@show inf_timing
# `total_time` is the total precompilation time across *all* loaded Julia modules
# `parcels` includes the *per-module* precompilation instructions and time 
total_time, parcels = SnoopCompile.parcel(inf_timing)
mktempdir() do tmpdir
    SnoopCompile.write(tmpdir, parcels)
    cp(joinpath(tmpdir, "precompile_BraketSimulator.jl"), joinpath(@__DIR__, "..", "src", "precompile.jl"); force=true)
    cp(joinpath(tmpdir, "precompile_BraketSimulatorPythonExt.jl"), joinpath(@__DIR__, "..", "ext", "BraketSimulatorPythonExt", "precompile.jl"); force=true)
end
