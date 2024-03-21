if ccall(:jl_generating_output, Cint, ()) == 1   # if we're precompiling the package
    let
        pybraket = pyimport("braket.circuits")
        n_qubits = 10
        c = pybraket.Circuit()
        c.h(0)
        # ghz
        for q in 1:n_qubits-1
            c.cnot(0, q)
        end
        c.state_vector()
        svs = StateVectorSimulator(0, 0)
        svs([c], shots=0)
    end
end