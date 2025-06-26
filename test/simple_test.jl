using BraketSimulator

# Try to access the functions
println("BranchedSimulator defined: ", isdefined(BraketSimulator, :BranchedSimulator))
println("calculate_probability_threshold defined: ", isdefined(BraketSimulator, :calculate_probability_threshold))
println("prune_paths! defined: ", isdefined(BraketSimulator, :prune_paths!))
println("allocate_shots defined: ", isdefined(BraketSimulator, :allocate_shots))
println("get_measurement_probabilities defined: ", isdefined(BraketSimulator, :get_measurement_probabilities))
println("apply_projection! defined: ", isdefined(BraketSimulator, :apply_projection!))
println("_apply_reset! defined: ", isdefined(BraketSimulator, :_apply_reset!))
println("measure_qubit! defined: ", isdefined(BraketSimulator, :measure_qubit!))
println("get_measurement_result defined: ", isdefined(BraketSimulator, :get_measurement_result))
println("set_variable! defined: ", isdefined(BraketSimulator, :set_variable!))
println("get_variable defined: ", isdefined(BraketSimulator, :get_variable))
println("new_to_circuit defined: ", isdefined(BraketSimulator, :new_to_circuit))
println("_evolve_branched_ast defined: ", isdefined(BraketSimulator, :_evolve_branched_ast))
