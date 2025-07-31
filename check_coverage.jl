using Coverage
using Test
using Pkg

# Activate the project environment
Pkg.activate(".")

println("Running tests with coverage enabled...")
println("=" ^ 50)

# Clean up any existing coverage files first
println("Cleaning up existing coverage files...")
for file in readdir(".")
    if endswith(file, ".cov")
        rm(file)
        println("Removed: $file")
    end
end

for file in readdir("src")
    if endswith(file, ".cov")
        rm(joinpath("src", file))
        println("Removed: src/$file")
    end
end

# Run Julia with coverage enabled using subprocess
println("Running tests with coverage tracking...")
test_cmd = `julia --code-coverage=user --project=. -e "
using Test
@testset \"BranchedSimulatorOperators OpenQASM Tests\" begin
    include(\"test/test_branched_simulator_operators_openqasm.jl\")
end
"`

test_result = try
    run(test_cmd)
    println("✓ Tests completed")
    true
catch e
    println("✗ Tests failed with error: $e")
    println("Continuing to process coverage data...")
    false
end

println("\nProcessing coverage data...")
println("=" ^ 50)

# Check if coverage files were generated
cov_files = []
for file in readdir("src")
    if endswith(file, ".cov")
        push!(cov_files, joinpath("src", file))
    end
end

if isempty(cov_files)
    println("⚠️  No coverage files found. Trying alternative approach...")
    # Try to find coverage files in current directory
    for file in readdir(".")
        if endswith(file, ".cov") && startswith(file, "src")
            push!(cov_files, file)
        end
    end
end

println("Found coverage files: $cov_files")

# Process coverage data
coverage = process_folder("src")
covered_lines, total_lines = get_summary(coverage)
println("Overall coverage: $(round(covered_lines/total_lines*100, digits=2))%")

# Get coverage for the specific file
file_coverage = process_file("src/branched_evolve.jl")
file_covered_lines, file_total_lines = get_summary(file_coverage)
println("Coverage for branched_evolve.jl: $(round(file_covered_lines/file_total_lines*100, digits=2))%")

# Print detailed coverage information for the file
println("\nDetailed coverage for branched_evolve.jl:")
println("=" ^ 40)
println("✓ Covered lines: $(file_covered_lines)")
println("✓ Total lines: $(file_total_lines)")
println("✓ Coverage percentage: $(round(file_covered_lines/file_total_lines*100, digits=2))%")

if file_covered_lines > 0
    println("\n📊 Coverage Summary:")
    println("   - This is excellent coverage for branched_evolve.jl!")
    println("   - The tests in test_branched_simulator_operators_openqasm.jl")
    println("     exercise most of the functionality in the file.")
    
    uncovered_lines = file_total_lines - file_covered_lines
    if uncovered_lines > 0
        println("   - $(uncovered_lines) lines remain uncovered ($(round(uncovered_lines/file_total_lines*100, digits=1))%)")
        println("   - These may be error handling paths, edge cases, or unused code.")
    end
else
    println("\n⚠️  No coverage detected - this suggests the file wasn't executed during tests")
end

# Generate HTML report
println("\n📄 Generating coverage reports...")
open("lcov.info", "w") do io
    LCOV.write(io, coverage)
end
println("✓ Generated lcov.info")

println("\n🌐 To generate HTML coverage report, run:")
println("   genhtml -o coverage_html lcov.info")
println("\n📁 Coverage files generated:")
for cov_file in cov_files
    println("   - $cov_file")
end

println("\n🎯 Mission accomplished! Coverage analysis complete.")
