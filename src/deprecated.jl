# All of the functions that take in a Program type (JAQCD) since it is no longer used

"""
	_prepare_program(circuit_ir::Program, inputs::Dict{String, <:Any}, shots::Int) -> (Program, Int)

Apply any `inputs` provided for the simulation. Return the `Program`
(with bound parameters) and the qubit count of the circuit.
"""
function _prepare_program(circuit_ir::Program, inputs::Dict{String, <:Any}, shots::Int) # nosemgrep
	operations::Vector{Instruction} = circuit_ir.instructions
	symbol_inputs = Dict(Symbol(k) => v for (k, v) in inputs)
	operations = [bind_value!(operation, symbol_inputs) for operation in operations]
	bound_program = Program(circuit_ir.braketSchemaHeader,
		operations,
		circuit_ir.results,
		circuit_ir.basis_rotation_instructions,
	)
	return bound_program, qubit_count(circuit_ir)
end

function _compute_results(simulator, program::Program, n_qubits, shots)
	analytic_results = shots == 0 && !isnothing(program.results) && !isempty(program.results)
	if analytic_results
		return _compute_exact_results(simulator, program, n_qubits)
	else
		return ResultTypeValue[]
	end
end
function _validate_circuit_ir(simulator, circuit_ir::Program, qubit_count::Int, shots::Int)
	_validate_ir_results_compatibility(simulator, circuit_ir.results, Val(:JAQCD))
	_validate_ir_instructions_compatibility(simulator, circuit_ir, Val(:JAQCD))
	_validate_shots_and_ir_results(shots, circuit_ir.results, qubit_count)
	return
end

function qubits(p::Program)
    inst_qubits = mapreduce(ix->ix.target, union, p.instructions, init=Set{Int}())
    bri_qubits  = mapreduce(ix->ix.target, union, p.basis_rotation_instructions, init=Set{Int}())
    res_qubits  = mapreduce(ix->(hasproperty(ix, :targets) && !isnothing(ix.targets)) ? reduce(vcat, ix.targets) : Set{Int}(), union, p.results, init=Set{Int}())
    return union(inst_qubits, bri_qubits, res_qubits)
end

qubit_count(p::Program) = length(qubits(p))

function Base.convert(::Type{Program}, c::Circuit) # nosemgrep
    lowered_rts = map(StructTypes.lower, c.result_types)
    header = braketSchemaHeader("braket.ir.jaqcd.program" ,"1")
    return Program(header, c.instructions, lowered_rts, c.basis_rotation_instructions)
end
Program(c::Circuit) = convert(Program, c)

