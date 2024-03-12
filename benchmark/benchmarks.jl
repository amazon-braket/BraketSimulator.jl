using BraketSimulator, Braket, PythonCall, BenchmarkTools

using Braket: Instruction, NoiseModel, add_noise!, GateCriteria

nm = NoiseModel() 
add_noise!(nm, BitFlip(0.1), GateCriteria(H))
add_noise!(nm, PhaseFlip(0.1), GateCriteria(CNot))

max_qubits(::Val{:noise}) = 16
max_qubits(::Val{:pure})  = 32

gate_operations  = pyimport("braket.default_simulator.gate_operations")
noise_operations = pyimport("braket.default_simulator.noise_operations")
local_sv         = pyimport("braket.default_simulator.state_vector_simulation")
local_dm         = pyimport("braket.default_simulator.density_matrix_simulation")
qml              = pyimport("pennylane")
np               = pyimport("numpy")
pnp              = pyimport("pennylane.numpy")
nx               = pyimport("networkx")

suite = BenchmarkGroup()
include("gate_kernels.jl")
include("qaoa.jl")
include("qft.jl")
include("ghz.jl")
#include("vqe.jl")

# this is expensive! only do it if we're sure we need to regen parameters
if !isfile("params.json")
    tune!(suite)
    BenchmarkTools.save("params.json", params(suite))
end
loadparams!(suite, BenchmarkTools.load("params.json")[1], :evals, :samples);
results = run(suite, verbose = true)
BenchmarkTools.save("results.json", results)
