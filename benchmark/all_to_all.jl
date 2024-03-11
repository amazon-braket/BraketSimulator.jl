using Braket, BraketStateVector

qubit_count = 22
n_layers = qubit_count^2
instructions = Braket.Instruction[]
for layer = 1:n_layers
    if isodd(layer)
        # CNots
        for qubit = 0:2:qubit_count-2
            push!(instructions, Braket.Instruction(CNot(), [qubit, qubit + 1]))
        end
    else
        # Hadamards 
        for qubit = 0:qubit_count-1
            push!(instructions, Braket.Instruction(H(), [qubit]))
        end
    end
end
simulation = StateVectorSimulator(qubit_count, 0)
simulation = evolve!(simulation, instructions)

simulation = StateVectorSimulator(qubit_count, 0)
p = @timev begin
    simulation = evolve!(simulation, instructions)
    nothing
end
