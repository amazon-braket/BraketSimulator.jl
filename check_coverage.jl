using Coverage
using Test
using Pkg

# Activate the test environment
Pkg.activate("test")

# Run only the specific test file with coverage enabled
try
    @testset "BranchedSimulatorOperators OpenQASM Tests" begin
        include("test/test_branched_simulator_operators_openqasm.jl")
        include("test/test_aliasing.jl")
    end
catch e
    println("Tests failed with error: $e")
    println("Continuing to process coverage data...")
end

# Process coverage data
coverage = process_folder("src")
covered_lines, total_lines = get_summary(coverage)
println("Overall coverage: $(round(covered_lines/total_lines*100, digits=2))%")

# Get coverage for the specific file
file_coverage = process_file("src/branched_evolve_operators.jl")
file_covered_lines, file_total_lines = get_summary(file_coverage)
println("Coverage for branched_evolve_operators.jl: $(round(file_covered_lines/file_total_lines*100, digits=2))%")

# Print detailed coverage information for the file
println("\nDetailed coverage for branched_evolve_operators.jl:")
# Skip detailed line-by-line coverage for now as it's causing issues
println("File has $(file_covered_lines) covered lines out of $(file_total_lines) total lines")

# Generate HTML report
open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end

println("\nTo generate HTML coverage report, run:")
println("genhtml -o coverage_html lcov.info")
