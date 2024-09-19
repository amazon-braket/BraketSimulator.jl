const ALL_QUBITS = -1

"""
    Circuit

A representation of a quantum circuit that contains the instructions to be performed on a
quantum device and the requested result types.

See:
  - [Gates](@ref) for all of the supported gates. 
  - [Noises](@ref) for all of the supported noise operations.
  - [Results](@ref) for all of the supported result types.
"""
mutable struct Circuit
    instructions::Vector{Instruction{<:Operator}}
    result_types::Vector{Result}
    basis_rotation_instructions::Vector{Instruction{<:Operator}}
    qubit_observable_mapping::Dict{Int, Observables.Observable}
    qubit_observable_target_mapping::Dict{Int, Tuple}
    qubit_observable_set::Set{Int}
    parameters::Set{FreeParameter}
    observables_simultaneously_measureable::Bool
    measure_targets::Vector{Int} 

    @doc """
        Circuit()
    
    Construct an empty `Circuit`.
    """
    Circuit() = new([], [], [], Dict(), Dict(), Set{Int}(), Set{FreeParameter}(), true, Int[])
end

"""
    qubits(c::Circuit) -> QubitSet

Returns a [`QubitSet`](@ref) containing all qubits
that `c` is defined on.

# Examples
```jldoctest
julia> c = Circuit();

julia> add_instruction!(c, Instruction(H(), 0));

julia> add_instruction!(c, Instruction(CNot(), [0, 1]));

julia> qubits(c)
QubitSet with 2 elements:
  0
  1
```
"""
qubits(c::Circuit) = (qs = union!(mapreduce(ix->ix.target, union, c.instructions, init=Set{Int}()), c.qubit_observable_set); QubitSet(qs))
function qubits(p::Program)
    inst_qubits = mapreduce(ix->ix.target, union, p.instructions, init=Set{Int}())
    bri_qubits  = mapreduce(ix->ix.target, union, p.basis_rotation_instructions, init=Set{Int}())
    res_qubits  = mapreduce(ix->(hasproperty(ix, :targets) && !isnothing(ix.targets)) ? reduce(vcat, ix.targets) : Set{Int}(), union, p.results, init=Set{Int}())
    return union(inst_qubits, bri_qubits, res_qubits)
end
"""
    qubit_count(c::Circuit) -> Int

Returns the number of qubits that `c` is defined on.

# Examples
```jldoctest
julia> c = Circuit();

julia> add_instruction!(c, Instruction(H(), 0));

julia> add_instruction!(c, Instruction(CNot(), [0, 1]));

julia> qubit_count(c)
2
```
"""
qubit_count(c::Circuit) = length(qubits(c))
qubit_count(p::Program) = length(qubits(p))

function Base.convert(::Type{Program}, c::Circuit) # nosemgrep
    lowered_rts = map(StructTypes.lower, c.result_types)
    header = braketSchemaHeader("braket.ir.jaqcd.program" ,"1")
    return Program(header, c.instructions, lowered_rts, c.basis_rotation_instructions)
end
Program(c::Circuit) = convert(Program, c)

extract_observable(rt::ObservableResult) = rt.observable
extract_observable(p::Probability) = Observables.Z()
extract_observable(rt::Result) = nothing 

function _encounter_noncommuting_observable!(c::Circuit)
    c.observables_simultaneously_measureable = false
    empty!(c.qubit_observable_mapping)
    empty!(c.qubit_observable_target_mapping)
    return c
end

function tensor_product_index_dict(o::Observables.TensorProduct, obs_targets::QubitSet)
    factors  = copy(o.factors)
    total    = qubit_count(first(factors))
    obj_dict = Dict{Int, Any}()
    i        = 0
    while length(factors) > 0
        if i >= total
            popfirst!(factors)
            if !isempty(factors)
                total += qubit_count(first(factors))
            end
        end
        if !isempty(factors)
            front = total - qubit_count(first(factors))
            obj_dict[i] = (first(factors), tuple([obs_targets[ii] for ii in front+1:total]...))
        end
        i += 1
    end
    return obj_dict
end

basis_rotation_gates(o::Observables.H) = (Ry(-π/4),)
basis_rotation_gates(o::Observables.X) = (H(),)
basis_rotation_gates(o::Observables.I) = ()
basis_rotation_gates(o::Observables.Z) = ()
basis_rotation_gates(o::Observables.Y) = (Z(), S(), H())
basis_rotation_gates(o::Observables.TensorProduct) = tuple(reduce(vcat, basis_rotation_gates.(o.factors))...)
basis_rotation_gates(o::Observables.HermitianObservable) = (Unitary(Matrix(adjoint(eigvecs(o.matrix)))),)

function fix_endianness(mat::Matrix)
    size(mat) != (4, 4) && return mat
    new_mat = copy(mat)
    for row in axes(new_mat, 1)
        idata  = new_mat[row, 2]
        new_mat[row, 2] = new_mat[row, 3]
        new_mat[row, 3] = idata
    end
    for col in axes(new_mat, 2)
        idata  = new_mat[2, col]
        new_mat[2, col] = new_mat[3, col]
        new_mat[3, col] = idata
    end
    return new_mat
end

function _observable_to_instruction(observable::Observables.Observable, target_list)::Vector{Instruction{<:Operator}}
    rotation_gates = collect(basis_rotation_gates(observable))
    return map(rotation_gates) do gate
        if gate isa Unitary && length(target_list) == 2
            return Instruction(Unitary(fix_endianness(gate.matrix)), target_list)
        else
            return Instruction(gate, target_list)
        end
    end
end

"""
    basis_rotation_instructions!(c::Circuit)

Gets a list of basis rotation instructions and stores them in the circuit `c`.
These basis rotation instructions are added if result types are requested for
an observable other than Pauli-Z.

This only makes sense if all observables are simultaneously measurable;
if not, this method will return an empty list.
"""
function basis_rotation_instructions!(c::Circuit)
    basis_rotation_instructions = Instruction[]
    all_qubit_observable = get(c.qubit_observable_mapping, ALL_QUBITS, nothing)
    if !isnothing(all_qubit_observable)
        c.basis_rotation_instructions = reduce(vcat, _observable_to_instruction(all_qubit_observable, target) for target in qubits(c))
        return c
    end
    unsorted = collect(Set(values(c.qubit_observable_target_mapping)))
    target_lists = sort(unsorted)
    for target_list in target_lists
        observable = c.qubit_observable_mapping[first(target_list)]
        append!(basis_rotation_instructions, _observable_to_instruction(observable, target_list))
    end
    c.basis_rotation_instructions = basis_rotation_instructions
    return c
end

function add_to_qubit_observable_mapping!(c::Circuit, o::Observables.Observable, obs_targets::QubitSet)
    targets = length(obs_targets) > 0 ? obs_targets : collect(c.qubit_observable_set)
    all_qubits_observable = get(c.qubit_observable_mapping, ALL_QUBITS, nothing)
    id = Observables.I()
    for ii in 1:length(targets)
        target             = targets[ii]
        new_observable     = o
        current_observable = !isnothing(all_qubits_observable) ? all_qubits_observable : get(c.qubit_observable_mapping, target, nothing)
        add_observable     = isnothing(current_observable) || (current_observable == id && new_observable != id)
        !add_observable && current_observable != id && new_observable != id && new_observable != current_observable && return _encounter_noncommuting_observable!(c)
        if !isempty(obs_targets)
            new_targets = tuple(obs_targets...)
            if add_observable
                c.qubit_observable_target_mapping[target] = new_targets
                c.qubit_observable_mapping[target] = new_observable
            elseif qubit_count(new_observable) > 1
                current_target = get(c.qubit_observable_target_mapping, target, nothing)
                !isnothing(current_target) && current_target != new_targets && _encounter_noncommuting_observable!(c)
            end
        end
    end
    if isempty(obs_targets)
        !isnothing(all_qubits_observable) && all_qubits_observable != o && return _encounter_noncommuting_observable!(c)
        c.qubit_observable_mapping[ALL_QUBITS] = o
    end
    return c
end

function add_to_qubit_observable_mapping!(c::Circuit, o::Observables.TensorProduct, obs_targets::QubitSet)
    targets = length(obs_targets) != 0 ? obs_targets : collect(c.qubit_observable_set)
    all_qubits_observable = get(c.qubit_observable_mapping, ALL_QUBITS, nothing)
    tensor_product_dict = length(targets) > 0 ? tensor_product_index_dict(o, QubitSet(targets)) : Dict()
    id = Observables.I()
    for ii in 1:length(targets)
        target             = targets[ii]
        new_observable     = tensor_product_dict[ii-1][1]
        current_observable = !isnothing(all_qubits_observable) ? all_qubits_observable : get(c.qubit_observable_mapping, target, nothing)
        add_observable     = isnothing(current_observable) || (current_observable == id && new_observable != id)
        !add_observable && current_observable != id && new_observable != id && new_observable != current_observable && return _encounter_noncommuting_observable!(c)
        if !isempty(obs_targets)
            new_targets = tensor_product_dict[ii-1][2]
            if add_observable
                c.qubit_observable_target_mapping[target] = new_targets
                c.qubit_observable_mapping[target] = new_observable
            elseif qubit_count(new_observable) > 1
                current_target = get(c.qubit_observable_target_mapping, target, nothing)
                !isnothing(current_target) && current_target != new_targets && return _encounter_noncommuting_observable!(c)
            end
        end
    end
    return c
end
add_to_qubit_observable_set!(c::Circuit, rt::ObservableResult) = union!(c.qubit_observable_set, Set(rt.targets))
# exclude AdjointGradient from coverage for now
# as we don't yet implement this, so don't have a test for it
# COV_EXCL_START
add_to_qubit_observable_set!(c::Circuit, rt::AdjointGradient)  = union!(c.qubit_observable_set, Set(reduce(union, rt.targets)))
# COV_EXCL_STOP
add_to_qubit_observable_set!(c::Circuit, rt::Result) = c.qubit_observable_set

function _check_if_qubit_measured(c::Circuit, qubit::Int)
    isempty(c.measure_targets) && return
    # check if there is a measure instruction on the targeted qubit(s)
    isempty(intersect(c.measure_targets, qubit)) || error("cannot apply instruction to measured qubits.")
end
_check_if_qubit_measured(c::Circuit, qubits) = foreach(q->_check_if_qubit_measured(c, Int(q)), qubits)

function add_instruction!(c::Circuit, ix::Instruction{O}) where {O<:Operator}
    _check_if_qubit_measured(c, ix.target)
    to_add = [ix]
    if ix.operator isa QuantumOperator && Parametrizable(ix.operator) == Parametrized()
        for param in parameters(ix.operator)
            union!(c.parameters, (param,))
        end
    end
    if ix.operator isa Measure
        append!(c.measure_targets, ix.target)
    end
    push!(c.instructions, ix)
    return c
end

function to_circuit(v::Quasar.QasmProgramVisitor)
    c = Circuit()
    foreach(v.instructions) do ix
        sim_op = StructTypes.constructfrom(QuantumOperator, ix)
        op     = isempty(ix.controls) ? sim_op : Control(sim_op, tuple(map(c->getindex(c, 2), ix.controls)...))
        sim_ix = Instruction(op, ix.targets)
        add_instruction!(c, sim_ix)
    end
    for rt in v.results
        sim_rt = StructTypes.constructfrom(Result, rt)
        obs    = extract_observable(sim_rt)
        if !isnothing(obs) && c.observables_simultaneously_measureable && !(rt isa AdjointGradient)
            add_to_qubit_observable_mapping!(c, obs, sim_rt.targets)
        end
        add_to_qubit_observable_set!(c, sim_rt)
        push!(c.result_types, sim_rt)
    end
    return c
end

# semgrep rules can't handle this macro properly yet
function to_circuit(qasm_source::String, inputs)
    input_qasm = if endswith(qasm_source, ".qasm") && isfile(qasm_source)
        read(qasm_source, String)
    else
        qasm_source
    end
    endswith(input_qasm, "\n") || (input_qasm *= "\n")
    parsed  = parse_qasm(input_qasm)
    visitor = QasmProgramVisitor(inputs)
    visitor(parsed)
    return to_circuit(visitor) 
end
to_circuit(qasm_source::String) = to_circuit(qasm_source, Dict{String, Float64}())
