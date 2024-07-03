using Coverage
# process '*.cov' files
coverage = process_folder()
coverage = append!(coverage, process_folder("ext/BraketSimulatorPythonExt/"))
coverage = append!(coverage, process_folder("ext/BraketSimulatorBraketExt/"))
coverage = merge_coverage_counts(coverage)
# Get total coverage for all Julia files
covered_lines, total_lines = get_summary(coverage)
@show covered_lines/total_lines

open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end

for fi in readdir("src")
    if endswith(fi, ".jl")
        println("Coverage for file $fi")
        @show get_summary(process_file(joinpath("src", fi)))
    end
end

# uncomment this if you have `genhtml` installed
# to generate HTML coverage info
#run(`genhtml -o .coverage/ lcov.info`)

clean_folder(".")
