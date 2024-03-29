function _check_annotations(node::OpenQASM3.Statement)
    if length(node.annotations) > 0
        @warn "Unsupported annotations $(node.annotations) at $(node.span)."
    end
    return nothing
end
_check_annotations(node) = nothing

value(exp::OpenQASM3.NothingExpression) = 1
value(exp) = exp.value

struct WalkerOutput
    ixs::Vector
    results::Vector{Braket.Result}
    WalkerOutput() = new([], Braket.AbstractProgramResult[])
end

abstract type AbstractQASMContext end
struct NothingContext <: AbstractQASMContext end
const ExternDef = OpenQASM3.ExternDeclaration

# We use `ConstOnlyVis` to prevent subroutines and defcals seeing non-const global QASM variables
@enum OpenQASMGlobalVisibility DefaultVis ConstOnlyVis
abstract type ClassicalVarKind end
struct ClassicalVar <: ClassicalVarKind end
struct ClassicalConst <: ClassicalVarKind end
struct ClassicalInput <: ClassicalVarKind end

mutable struct ClassicalDef{T, S, K}
    value::T
    valtype::S
    kind::K
end

mutable struct SubroutineDef{T, S, F, C}
    args::Vector{T}
    body_gen::F
    return_type::S
    ctx::C
end

mutable struct GateDef{T, Q, F, C}
    args::Vector{T}
    qubits::Vector{Q}
    body_gen::F
    ctx::C
end
Braket.qubit_count(gd::GateDef) = length(gd.qubits)
function (gd::GateDef)(arg_vals...)
    for (arg, val) in zip(gd.args, arg_vals)
        gd.ctx.definitions[arg] = ClassicalDef(val, typeof(val), ClassicalVar())
    end
    return gd.body_gen()
end
Base.show(io::IO, gd::GateDef) = print(io, "GateDef on qubits $(gd.qubits)")

struct QubitDef
    size::Int
end
QubitDef(::Nothing) = QubitDef(1)
const QASMDefinition = Union{ClassicalDef, GateDef, QubitDef, ExternDef, SubroutineDef}

Base.length(def::QubitDef) = def.size
Base.length(def::QASMDefinition) = length(def.value)

"""
    QASMGlobalContext{W}(;[ ext_lookup,])

The top-level context, which hold variables that live in the global QASM scope.
In addition to global variable definitions, this structure holds externally
defined input scalars, waveforms, and ports. These can be made accessible
to the OpenQASM program via `input` or `extern` statements.

Finally, there is a `visiblity` field that determines whether non-const
definitions should be hidden or visible. This is used to hide non-const global
variables when processing subroutines, which can only refer to
outer-scope variables if they are global constants.

See also [`QASMBlockContext`](@ref).
"""
mutable struct QASMGlobalContext{W, F} <: AbstractQASMContext
    definitions::Dict{String, QASMDefinition}
    qubit_mapping::Dict{Union{String, Tuple{String, Int}}, Vector{Int}}
    n_qubits::Int
    ext_lookup::F
    visiblity::OpenQASMGlobalVisibility
    output::WalkerOutput
end

function QASMGlobalContext{W}(ext_lookup::F=Dict{String,Any}()) where {W, F}
    ctx = QASMGlobalContext{W, F}(
        Dict{String, QASMDefinition}(),
        Dict{Union{String, Tuple{String, Int}}, Vector{Int}}(),
        0,
        ext_lookup,
        DefaultVis,
        WalkerOutput()
    )
    return ctx
end

"""
    QASMBlockContext(parent)

Data for tracking variables in a local QASM block scope, such as inside
if-else or loop blocks.

See also [`QASMGlobalContext`](@ref).
"""
struct QASMBlockContext{P <: AbstractQASMContext} <: AbstractQASMContext
    parent::P
    qubit_mapping::Dict{Union{String, Tuple{String, Int}}, Vector{Int}}
    definitions::Dict{String, QASMDefinition}
    function QASMBlockContext(parent::P) where {P <: AbstractQASMContext}
        return new{P}(parent, parent.qubit_mapping, Dict{String, QASMDefinition}())
    end
end

"""
    QASMSubroutineContext(parent)

Data for tracking variables in a QASM subroutine scope.

See also [`QASMGlobalContext`](@ref).
"""
mutable struct QASMSubroutineContext{P <: AbstractQASMContext} <: AbstractQASMContext
    parent::P
    qubit_mapping::Dict{Union{String, Tuple{String, Int}}, Vector{Int}}
    definitions::Dict{String, QASMDefinition}
    function QASMSubroutineContext(parent::P) where {P <: AbstractQASMContext}
        return new{P}(parent, Dict{String, Int}(), Dict{String, QASMDefinition}())
    end
end

"""
    QASMGateDefContext(parent)

Data for tracking variables in a QASM gate definition scope.

See also [`QASMGlobalContext`](@ref).
"""
mutable struct QASMGateDefContext{P <: AbstractQASMContext} <: AbstractQASMContext
    parent::P
    qubit_mapping::Dict{Union{String, Tuple{String, Int}}, Vector{Int}}
    n_qubits::Int
    definitions::Dict{String, QASMDefinition}
    output::WalkerOutput
    function QASMGateDefContext(parent::P) where {P <: AbstractQASMContext}
        return new{P}(parent, Dict{String, Int}(), 0, Dict{String, QASMDefinition}(), WalkerOutput())
    end
end

Base.push!(ctx::AbstractQASMContext, ix::Braket.Instruction{O}) where {O<:Braket.Operator} = push!(output(ctx).ixs, ix)
Base.push!(ctx::AbstractQASMContext, rt::Braket.Result)      = push!(output(ctx).results, rt)

output(ctx::QASMGlobalContext)     = ctx.output
output(ctx::QASMGateDefContext)    = ctx.output
output(ctx::AbstractQASMContext)   = output(parent(ctx))

Base.parent(ctx::AbstractQASMContext) = ctx.parent
Base.parent(::QASMGlobalContext)      = NothingContext() 

const NumberLiteral = Union{OpenQASM3.IntegerLiteral, OpenQASM3.FloatLiteral, OpenQASM3.BooleanLiteral, OpenQASM3.BitstringLiteral}

Base.convert(::Type{Int},     v::OpenQASM3.IntegerLiteral) = v.value
Base.convert(::Type{T},       v::OpenQASM3.IntegerLiteral) where {T<:Real} = convert(T, v.value)
Base.convert(::Type{Float64}, v::OpenQASM3.FloatLiteral) = Float64(v.value)
Base.convert(::Type{Float32}, v::OpenQASM3.FloatLiteral) = Float32(v.value)
Base.convert(::Type{Float16}, v::OpenQASM3.FloatLiteral) = Float16(v.value)

is_defined(name::String, ctx::AbstractQASMContext; local_only=false) = return haskey(ctx.definitions, name) || (!local_only && is_defined(name, parent(ctx)))
is_defined(name::String, ::Nothing; kwargs...) = false
is_defined(name::String, ::NothingContext; kwargs...) = false
is_defined(name, ctx::AbstractQASMContext; kwargs...) = false

function interpret!(output::WalkerOutput, node::T) where {T}
    @error "Unsupported instruction at $(node.span): $T"
    return nothing
end

_check_type(::Type{T}, val::T, name::String, span) where {T} = val
_check_type(::Type{T}, val::V, name::String, span) where {T, V} = error("Expected `$name` to be type $T but got $V, at $span.")

"""
    lookup_def(::Type, name, context::AbstractQASMContext; span)

Lookup a definition by recursing through the hierarchy of contexts. If an
identifier is not defined in the context for the current scope, try the parent.
"""
function lookup_def(::Type{T}, name::String, def::D, ctx::QASMGlobalContext; span=()) where {T, D<:QASMDefinition}
    ctx.visiblity == ConstOnlyVis && !_is_global_const(def) && error("Attempt to use non-const global `$name` in subroutine or defcal.")
    return _check_type(T, def, name, span)
end
lookup_def(::Type{T}, name::String, def::D, ctx::AbstractQASMContext; span=()) where {T, D<:QASMDefinition} = return _check_type(T, def, name, span)
function lookup_def(::Type{GateDef}, name::String, ctx::AbstractQASMContext; span=())::GateDef
    haskey(ctx.definitions, name) && return _check_type(GateDef, ctx.definitions[name], name, span)
    return lookup_def(GateDef, name, parent(ctx); span=span)
end
function lookup_def(::Type{T}, name::String, ctx::AbstractQASMContext; span=()) where {T}
    def = haskey(ctx.definitions, name) ? ctx.definitions[name] : lookup_def(T, name, parent(ctx); span=span)
    lookup_def(T, name, def, ctx, span=span)
end
lookup_def(::Type{T}, name::String, def::Nothing, ctx::AbstractQASMContext; span=()) where {T} = throw("definition $name of type $T not found in global context!")
lookup_def(::Type{T}, name::String, ctx::NothingContext; span=()) where {T} = throw("definition $name of type $T not found in global context!")

id(node::OpenQASM3.Identifier)            = node.name
id(node::OpenQASM3.QASMNode)              = id(node.name)
id(node::OpenQASM3.IndexExpression)       = id(node.collection)
id(node::OpenQASM3.IODeclaration)         = id(node.identifier)
id(node::OpenQASM3.ExternDeclaration)     = id(node.identifier)
id(node::OpenQASM3.ClassicalDeclaration)  = id(node.identifier)
id(node::OpenQASM3.ClassicalAssignment)   = id(node.lvalue)
id(node::OpenQASM3.QubitDeclaration)      = id(node.qubit)
id(node::OpenQASM3.ConstantDeclaration)   = id(node.identifier)
id(node::OpenQASM3.ForInLoop)             = id(node.identifier)
id(node::OpenQASM3.QuantumGateModifier)   = node.modifier
id(node::OpenQASM3.UnaryOperator{O}) where {O} = O
id(node::OpenQASM3.BinaryOperator{O}) where {O} = O
id(node::Tuple{String, <:Any})            = node[1]
id(node::String)                          = node

(ctx::AbstractQASMContext)(node::T) where {T<:NumberLiteral} = node.value
(ctx::AbstractQASMContext)(node::T) where {T<:Number}        = node
(ctx::AbstractQASMContext)(node::BitVector)                  = node
_check_def_value(def::Nothing, name::String, span)           = error("'$name' referenced in $span, but not defined.")
_check_def_value(def, name::String, span)                    = isnothing(def.value) ? error("Uninitialized variable `$name` used in $span.") : nothing

(ctx::AbstractQASMContext)(::Nothing) = nothing
function (ctx::AbstractQASMContext)(node::OpenQASM3.BitstringLiteral)
    val = falses(node.width)
    for ix in 0:node.width-1
        val[end - ix] = 1 & (node.value >> ix)
    end
    return val
end
function (ctx::AbstractQASMContext)(node::OpenQASM3.SizeOf)
    target = ctx(node.target)
    dim    = ctx(node.index)
    d      = isnothing(dim) ? 1 : dim
    return size(target, d)
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.Identifier)
    name = id(node)
    name == "pi"    || name == "π" && return π
    name == "euler" || name == "ℇ" && return ℯ 
    def  = lookup_def(QASMDefinition, name, ctx; span=node.span)
    if def isa ClassicalDef
        _check_def_value(def, name, node.span)
        # TODO: Maybe return a wrapper to track OpenQASM types, const status, etc.
        return ctx(def.value)
    elseif def isa QubitDef
        return resolve_qubit(name, ctx)
    end
end
(ctx::AbstractQASMContext)(node::OpenQASM3.Identifier, len::Int) = ctx(node)
(ctx::AbstractQASMContext)(node::OpenQASM3.ArrayLiteral) = [ctx(v) for v in node.values]
(ctx::AbstractQASMContext)(node::OpenQASM3.DiscreteSet)  = [ctx(v) for v in node.values]

function (ctx::AbstractQASMContext)(node::OpenQASM3.RangeDefinition)
    start = ctx(node.start)
    stop  = ctx(node.stop)
    step  = isnothing(node.step) ? 1 : ctx(node.step)
    return range(start, step=step, stop=stop)
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.RangeDefinition, len::Int)
    start = ctx(node.start)
    stop  = ctx(node.stop)
    if !isnothing(stop)
        stop < 0 && (stop = len + stop)
    else
        stop = len-1
    end
    step  = isnothing(node.step) ? 1 : ctx(node.step)
    return range(start, step=step, stop=stop)
end

_get_indices(node::OpenQASM3.RangeDefinition{R, S, T}, name::AbstractString, ctx::AbstractQASMContext) where {R, S, T} = ctx(node, length(lookup_def(QASMDefinition, name, ctx)))
_get_indices(node::Vector{OpenQASM3.RangeDefinition{R, S, T}}, name::AbstractString, ctx::AbstractQASMContext) where {R, S, T} = _get_indices(node[1], name, ctx)
_get_indices(node::T, name::AbstractString, ctx::AbstractQASMContext) where {T} = ctx(node)
function _lookup_name(node::OpenQASM3.IndexExpression, ctx::AbstractQASMContext)
    name    = id(node)
    indices = _get_indices(node.index, name, ctx)
    return (name, [index for index in indices])
end
_lookup_name(node, ctx) = id(node)

_index_val(val, inds, typ) = val[inds .+ 1]
function _index_val(val::T, inds, typ) where {T<:Integer}
    width = sizeof(T)
    v_digits = digits(val, base=2, pad=OpenQASM3.value(typ.size))
    vals = reverse(v_digits)[inds .+ 1]
    return vals
end
_generate_elements_from_obj(obj::ClassicalDef, name::String, inds) = _index_val(obj.value, inds, obj.valtype)
_generate_elements_from_obj(obj::QubitDef, name::String, inds)       = [(name, inds)]
function (ctx::AbstractQASMContext)(node::OpenQASM3.IndexExpression)
    names    = _lookup_name(node, ctx)
    obj_name = names[1]
    mapped_pairs = map(names[2]) do obj_inds
        obj = lookup_def(QASMDefinition, obj_name, ctx; span=node.span)
        return _generate_elements_from_obj(obj, obj_name, obj_inds)
    end
    return collect(Iterators.flatten(mapped_pairs))
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.ExpressionStatement)
    _check_annotations(node)
    ctx(node.expression)
    return
end

(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:+}}) = ctx(node.lhs) + ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:-}}) = ctx(node.lhs) - ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:*}}) = ctx(node.lhs) * ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:/}}) = ctx(node.lhs) / ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:&}}) = (&).(ctx(node.lhs), ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:|}}) = (|).(ctx(node.lhs), ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:^}}) = (⊻).(ctx(node.lhs), ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:<<}}) = ctx(node.lhs) << ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:>>}}) = ctx(node.lhs) >> ctx(node.rhs)
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{Symbol("**")}}) = ctx(node.lhs) ^ ctx(node.rhs)
    
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:<}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) < ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:>}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) > ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:(!=)}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) != ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:(==)}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) == ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:(<=)}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) <= ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{Val{:(>=)}}) = scalar_cast(OpenQASM3.BoolType, ctx(node.lhs) >= ctx(node.rhs))
(ctx::AbstractQASMContext)(node::OpenQASM3.BinaryExpression{O}) where O = error("Binary op $O not yet implemented.")

(ctx::AbstractQASMContext)(node::OpenQASM3.Cast{T}) where {T} = scalar_cast(T, ctx(node.argument))

(ctx::AbstractQASMContext)(node::OpenQASM3.UnaryExpression{Val{:-}}) = -ctx(node.expression)
(ctx::AbstractQASMContext)(node::OpenQASM3.UnaryExpression{Val{:!}}) = !scalar_cast(OpenQASM3.BoolType, ctx(node.expression))
(ctx::AbstractQASMContext)(node::OpenQASM3.UnaryExpression{Val{:~}}) = (!).(ctx(node.expression))
(ctx::AbstractQASMContext)(node::OpenQASM3.UnaryExpression{O}) where {O} = error("Unary op $(node.op) not yet implemented.")
(ctx::AbstractQASMContext)(nodes::Vector{T}) where {T} = map(ctx, nodes)
(ctx::AbstractQASMContext)(nodes::Vector{IntegerLiteral}, len::Int) = map(ctx, nodes)
(ctx::AbstractQASMContext)(nodes::Vector{T}, len::Int) where {T} = map(ctx, nodes, len)

# FIXME: These do not enforce widths yet.
scalar_matches_type(::Type{OpenQASM3.FloatType}, val::T, err_str::String) where {T<:Real}      = return
scalar_matches_type(::Type{OpenQASM3.ComplexType}, val::T, err_str::String) where {T<:Complex} = return
scalar_matches_type(::Type{OpenQASM3.IntType}, val::T, err_str::String) where {T<:Integer}     = return
scalar_matches_type(::Type{OpenQASM3.BitType}, val::T, err_str::String) where {T<:Integer}     = return
scalar_matches_type(::Type{OpenQASM3.BitType}, val::BitArray, err_str::String)                 = return
scalar_matches_type(::Type{OpenQASM3.BoolType}, val::Bool, err_str::String)                    = return
scalar_matches_type(::Type{OpenQASM3.UintType}, val::T, err_str::String) where {T<:Unsigned}   = return
scalar_matches_type(::Type{OpenQASM3.UintType}, val::T, err_str::String) where {T<:Integer}    = scalar_matches_type(OpenQASM3.UintType, convert(Unsigned, val), err_str)
# FIXME: OpenQASM3.jl should not be putting unsigned ints into Ints.
scalar_matches_type(::Type{O}, val::T, err_str::String) where {O<:OpenQASM3.UintType,T<:Integer}    = (val >= 0 && error(err_str * " $T"); return)
scalar_matches_type(::Type{T}, val::V, err_str::String) where {T, V} = error(err_str * " $V, $T")
scalar_matches_type(t::T, val::V, err_str::String) where {T, V} = scalar_matches_type(T, val, err_str * " $V")

scalar_matches_type(t::OpenQASM3.ArrayType, val::Vector{V}, err_str::String) where {V} = foreach(v->scalar_matches_type(t.base_type, v, err_str), val)

# FIXME: should check actual size
scalar_matches_type(::Type{OpenQASM3.FloatType}, val::OpenQASM3.FloatLiteral, err_str::String) = return
scalar_matches_type(::Type{OpenQASM3.IntType}, val::OpenQASM3.IntegerLiteral, err_str::String) = return

scalar_cast(::OpenQASM3.FloatType, val) = real(1.0 * val)
scalar_cast(::OpenQASM3.ComplexType, val) = complex(val)
scalar_cast(::OpenQASM3.IntType, val) = Int(val)
scalar_cast(::OpenQASM3.UintType, val) = UInt(val)
scalar_cast(::OpenQASM3.BoolType, val) = Bool(val)
scalar_cast(::OpenQASM3.BoolType, val::BitVector) = any(i->i>0, val)
scalar_cast(::Type{OpenQASM3.BoolType}, val) = Bool(val)
scalar_cast(::Type{OpenQASM3.BoolType}, val::Vector{Int}) = any(i->i>0, val)
scalar_cast(::Type{OpenQASM3.BoolType}, val::BitVector) = any(i->i>0, val)

_new_inner_scope(ctx::AbstractQASMContext)      = QASMBlockContext(ctx)
_new_subroutine_scope(ctx::AbstractQASMContext) = QASMSubroutineContext(ctx)
_new_gate_def_scope(ctx::AbstractQASMContext)   = QASMGateDefContext(ctx)

function get_pragma_arg(body::String)
    match_arg = match(r"\(([^\)]+)\)", body)
    stripped_body = replace(body, r"\(([^\)]+)\)"=>"")
    arg = isnothing(match_arg) ? nothing : match_arg.match 
    return stripped_body, arg
end
get_pragma_arg(body) = get_pragma_arg(String(body))
function (ctx::AbstractQASMContext)(::Val{:adjoint_gradient}, body::AbstractString)
    chopped_body = replace(chopprefix(body, "adjoint_gradient "), "\""=>"")
    segment_start = findfirst('(', chopped_body)
    segment_end   = findlast(')', chopped_body)
    # handle arg as OQ3 string 
    op_chunk = split(chopped_body[segment_start+1:segment_end-1], '+')
    ops = []
    targets = Vector{Int}[]
    for chunk in op_chunk
        op, qubit = ctx(Val(:operator), String(chunk))
        push!(ops, op) 
        push!(targets, qubit) 
    end
    params_str = chopped_body[segment_end+1:end]
    params     = filter(!isempty, String.(split(params_str, " ")))
    push!(ctx, Braket.AdjointGradient(sum(ops), targets, params))
    return
end

function (ctx::AbstractQASMContext)(::Val{:amplitude}, body::AbstractString)
    chopped_body = replace(chopprefix(body, "amplitude "), "\""=>"", "'"=>"")
    states       = split(chopped_body, " ")
    push!(ctx, Braket.Amplitude([String(replace(s, ","=>"")) for s in states]))
    return
end

function (ctx::AbstractQASMContext)(::Val{:state_vector}, body::AbstractString)
    push!(ctx, Braket.StateVector())
    return
end

function (ctx::AbstractQASMContext)(::Val{:density_matrix}, body::AbstractString)
    chopped_body = chopprefix(body, "density_matrix")
    targets = ctx(Val(:pragma), Val(:qubits), chopped_body)
    push!(ctx, Braket.DensityMatrix(targets))
    return
end

function (ctx::AbstractQASMContext)(::Val{:probability}, body::AbstractString)
    chopped_body = chopprefix(body, "probability")
    targets      = ctx(Val(:pragma), Val(:qubits), chopped_body)
    push!(ctx, Braket.Probability(targets))
    return
end

function parse_hermitian(arg)
    clean_arg = chop(arg, head=2, tail=2) # get rid of brackets on either side
    vecs      = split(replace(clean_arg, "["=>""), "],")
    h_mat     = Matrix{ComplexF64}(undef, length(vecs), length(vecs))
    for (i, v) in enumerate(vecs)
        for (j, elem) in enumerate(split(v, ","))
            h_mat[i,j] = tryparse(ComplexF64, elem)
        end
    end
    mat = Braket.Observables.HermitianObservable(h_mat)
    return mat
end

function parse_standard_op(str::AbstractString)
    if occursin('*', str)
        clean_str = replace(str, " "=>"")
        coeff = tryparse(Float64, split(clean_str, '*')[1])
        op = Braket.StructTypes.constructfrom(Braket.Observables.Observable, string(split(clean_str, '*')[end]))
    elseif occursin('-', str)
        coeff = -1.0
        op = Braket.StructTypes.constructfrom(Braket.Observables.Observable, string(split(str, '-')[end]))
    else
        coeff = 1.0
        op = Braket.StructTypes.constructfrom(Braket.Observables.Observable, string(str[1]))
    end
    return coeff * op
end

function parse_individual_op(op_str::AbstractString, ctx::AbstractQASMContext)
    head_chop = startswith(op_str, ' ') ? 1 : 0
    tail_chop = endswith(op_str, ' ') ? 1 : 0
    clean_op  = chop(op_str, head=head_chop, tail=tail_chop)
    stripped_op, arg = get_pragma_arg(clean_op)
    clean_arg        = isnothing(arg) ? nothing : replace(arg, "("=>"", ")"=>"")
    is_hermitian     = startswith(stripped_op, "hermitian")
    op = is_hermitian ? parse_hermitian(clean_arg) : parse_standard_op(stripped_op)
    qubits = nothing
    if is_hermitian || !isnothing(arg)
        qubit_str = is_hermitian ? replace(stripped_op, "hermitian "=>"") : clean_arg
        qubits    = ctx(Val(:pragma), Val(:qubits), qubit_str)
    end
    return op, qubits
end

function (ctx::AbstractQASMContext)(::Val{:operator}, op_str::String)
    is_tensor_prod = occursin('@', op_str)
    if is_tensor_prod
        op = Braket.Observables.Observable[]
        qubits = Int[]
        for op_ in split(op_str, '@')
            this_op, this_qubits = parse_individual_op(op_, ctx)
            push!(op, this_op)
            append!(qubits, this_qubits)
        end
        return Braket.Observables.TensorProduct(op), qubits
    else
        return parse_individual_op(op_str, ctx)
    end
end

lookup_mapping(q::String, ctx) = ctx.qubit_mapping[q]
lookup_mapping(q::Tuple{String, T}, ctx) where {T} = [ctx.qubit_mapping[q[1]][q_+1] for q_ in q[2]]
function (ctx::AbstractQASMContext)(::Val{:pragma}, ::Val{:qubits}, body::AbstractString)
    (occursin("all", body) || isempty(body)) && return nothing
    has_brakets = occursin('[', body)
    if has_brakets # more complicated...
        qubits = Int[]
        segment_start = 1
        segment_end   = findnext(']', body, segment_start)
        while !isnothing(segment_end)
            raw_chunk = chopprefix(String(body[segment_start:segment_end]), ",")
            oq3_chunk = raw_chunk * ";\n"
            parsed_chunk = OpenQASM3.parse(oq3_chunk)
            for node in parsed_chunk.statements
                expr = ctx(node.expression)
                append!(qubits, Iterators.flatten([lookup_mapping(q, ctx) for q in expr]))
            end
            segment_start = segment_end + 1
            !isnothing(segment_start) && (segment_start += 1)
            segment_end = !isnothing(segment_start) && segment_start <= length(body) ? findnext(']', body, segment_start) : nothing
        end
    else
        qubits = map(split(body, ",")) do q_name
            clean_name = String(replace(q_name, " "=>""))
            return collect(Iterators.flatten(ctx.qubit_mapping[q] for q in resolve_qubit(clean_name, ctx)))
        end
    end
    targets    = collect(Iterators.flatten(qubits))
    return targets
end

for (tag, type, type_str) in ((Val{:expectation}, :(Braket.Expectation), "expectation"), (Val{:variance}, :(Braket.Variance), "variance"), (Val{:sample}, :(Braket.Sample), "sample"))
    @eval begin
        function (ctx::AbstractQASMContext)(::$tag, body::AbstractString)
            chopped_body = String(chopprefix(body, $type_str * " "))
            op_str       = chopped_body
            op, targets  = ctx(Val(:operator), op_str)
            rt           = $type(op, targets)
            push!(ctx, rt)
            return
        end
    end
end

function (ctx::AbstractQASMContext)(::Val{:result}, body::String)
    chopped_body = chopprefix(body, "result ")
    result_type  = Val(Symbol(split(chopped_body, " ")[1]))
    ctx(result_type, chopped_body)
    return
end

function get_noise_arg(::Val{:kraus}, arg::AbstractString)
    stripped_arg = replace(arg, ")"=>"]", "("=>"[")
    raw_arg      = eval(Meta.parse(stripped_arg))
    noise_arg    = map(a->reduce(hcat, a), raw_arg)
    return noise_arg
end
function get_noise_arg(::Val{:noise}, arg::AbstractString)
    args = map(split(arg, ",")) do arg_str
        return tryparse(Float64, replace(arg_str, " "=>"", ")"=>"", "("=>""))
    end
    return args
end
function (ctx::AbstractQASMContext)(::Val{:noise}, body::String)
    stripped_body, arg = get_pragma_arg(String(chopprefix(body, "noise ")))
    split_body         = filter(!isempty, split(stripped_body, " "))
    raw_op             = split_body[1]
    rebuilt_qubits     = join(split_body[2:end], " ")
    targets            = ctx(Val(:pragma), Val(:qubits), rebuilt_qubits)
    noise_id           = Symbol(raw_op)
    noise_type         = Braket.StructTypes.subtypes(Noise)[noise_id]
    if raw_op == "kraus"
        noise_arg = get_noise_arg(Val(:kraus), arg)
        op        = Kraus(noise_arg)
        push!(ctx, Instruction{Kraus}(op, collect(Iterators.flatten(targets))))
    else
        args = get_noise_arg(Val(:noise), arg)
        op   = noise_type(args...)
        push!(ctx, Instruction{noise_type}(op, collect(Iterators.flatten(targets))))
    end
    return
end

function (ctx::AbstractQASMContext)(::Val{:unitary}, body::String)
    stripped_body, arg = get_pragma_arg(String(chopprefix(body, "unitary ")))
    split_body         = filter(!isempty, split(stripped_body, " "))
    rebuilt_qubits     = join(split_body[2:end], " ")
    targets            = ctx(Val(:pragma), Val(:qubits), rebuilt_qubits)
    cleaned_arg        = replace(arg, "("=>"", ")"=>"")
    raw_arg            = eval(Meta.parse(cleaned_arg))
    gate_arg           = reduce(hcat, raw_arg)
    op                 = Unitary(gate_arg)
    push!(ctx, Instruction{Unitary}(op, targets))
    return
end

function (ctx::AbstractQASMContext)(::Val{:pragma}, cmd::String)
    !startswith(cmd, "braket ") && error("pragma `$cmd` must begin with `braket `")
    occursin("verbatim", cmd) && return
    pragma_body = String(chopprefix(cmd, "braket "))
    end_char    = if occursin(' ', pragma_body) && occursin('(', pragma_body)
                      min(findfirst(' ', pragma_body)-1, findfirst('(', pragma_body)-1)
                  elseif occursin(' ', pragma_body) && !occursin('(', pragma_body)
                      findfirst(' ', pragma_body)-1
                  elseif !occursin(' ', pragma_body) && occursin('(', pragma_body)
                      findfirst('(', pragma_body)-1
                  end
    pragma_type = pragma_body[1:end_char]
    op_type     = Val(Symbol(pragma_type))
    ctx(op_type, pragma_body)
    return
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.Box)
    foreach(ctx, node.body)
    return nothing
end
(ctx::AbstractQASMContext)(node::OpenQASM3.Pragma) = ctx(Val(:pragma), node.command)

function _lookup_ext_scalar(name, ctx::QASMGlobalContext)
    isnothing(ctx.ext_lookup) && error("Context lookup is empty")
    val = get(ctx.ext_lookup, name, nothing)
    isnothing(val) || val isa Number || error("QASM program expects '$name' to be a scalar. Got '$val'.")
    return val
end
_lookup_ext_scalar(name, ctx::AbstractQASMContext) = _lookup_ext_scalar(name, parent(ctx))

# NOTE: `node.type` cannot be inferred. Could avoid runtime dispatch here
#       by branching on allowed types (i.e. manual dispatch).
function (ctx::AbstractQASMContext)(node::N, args...) where {N <: OpenQASM3.QASMNode}
    _check_annotations(node)
    name = id(node)
    _check_undefined(node, name, ctx)
    return ctx(node, name, args...)
end

(ctx::AbstractQASMContext)(node::OpenQASM3.NothingExpression) = nothing
(ctx::AbstractQASMContext)(node::OpenQASM3.IODeclaration{T, OpenQASM3.output}, name::String) where {T} = error("Output not supported.")
(ctx::AbstractQASMContext)(node::OpenQASM3.IODeclaration{T, IT}, name::String) where {T, IT} = error("IO type $IT not supported at $(node.span).")
function (ctx::AbstractQASMContext)(node::OpenQASM3.IODeclaration{T, OpenQASM3.input}, name::String) where {T<:OpenQASM3.ClassicalType}
    val = _lookup_ext_scalar(name, ctx)
    isnothing(val) && error("Input variable $name was not supplied at $(node.span).")
    scalar_matches_type(T, val, "Input variable $name at $(node.span): type does not match ")
    ctx.definitions[name] = ClassicalDef(val, node.type, ClassicalInput)
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.ExternDeclaration, name::String)
    ctx.definitions[name] = node
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.QuantumPhase)
    _check_annotations(node)
    if isempty(node.qubits)
        gate_qubits = collect(0:ctx.n_qubits-1)
    else
        gate_qubits = collect(Iterators.flatten(_get_qubits_from_mapping(node.qubits, ctx)))
    end
    gate_arg         = convert(Float64, ctx(node.argument))
    braket_gate_type = MultiQubitPhaseShift{length(gate_qubits)}
    braket_gate      = braket_gate_type(gate_arg...)
    modified_gate    = ctx(node.modifiers, braket_gate)
    push!(ctx, Instruction(modified_gate, gate_qubits))
    return
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.QuantumGateDefinition, name::String)
    gate_ctx = _new_gate_def_scope(ctx)
    args     = map(ctx, node.arguments)
    for (q_ix, q) in enumerate(node.qubits)
        q_name = id(q)
        gate_ctx.definitions[q_name] = QubitDef(1)
        gate_ctx.qubit_mapping[q_name] = [q_ix - 1]
        gate_ctx.qubit_mapping[(q_name, 0)] = [q_ix - 1]
        gate_ctx.n_qubits += 1
    end
    for statement in node.body
        gate_ctx(statement)
    end
    body_fn() = [Instruction(ix.operator, ix.target) for ix in gate_ctx.output.ixs]
    ctx.definitions[name] = GateDef(args, id.(node.qubits), body_fn, gate_ctx) 
    return
end

pad_qubits(op::Control, qc_diff::Int, target) = vcat(collect(0:qc_diff-1), target .+ qc_diff)
pad_qubits(op, qc_diff::Int, target) = target
n_arg(::Nothing) = 1
n_arg(arg)       = arg
(ctx::AbstractQASMContext)(node::QuantumGateModifier{OpenQASM3.Ctrl, E}, argument, n_arg::Int, g::G) where {G<:Braket.Gate, E<:OpenQASM3.Expression} = Control(g, ntuple(i->1, n_arg))
(ctx::AbstractQASMContext)(node::QuantumGateModifier{OpenQASM3.NegCtrl, E}, argument, n_arg::Int, g::G) where {G<:Braket.Gate, E<:OpenQASM3.Expression} = Control(g, ntuple(i->0, n_arg))
function (ctx::AbstractQASMContext)(mod::QuantumGateModifier{T, E}, argument, braket_gate::Vector{Instruction{O}}) where {T<:Union{OpenQASM3.Ctrl, OpenQASM3.NegCtrl}, O<:Braket.Operator, E<:OpenQASM3.Expression}
    braket_gate = map(braket_gate) do ix
        old_qc  = qubit_count(ix.operator)
        new_op  = ctx(mod, n_arg(argument), ix.operator)
        new_qc  = qubit_count(new_op)
        qc_diff = new_qc - old_qc
        new_ix  = Instruction(new_op, pad_qubits(new_op, qc_diff, ix.target))
    end
    return braket_gate
end
(ctx::AbstractQASMContext)(mod::QuantumGateModifier{T, E}, argument, braket_gate::G) where {T<:Union{OpenQASM3.Ctrl, OpenQASM3.NegCtrl}, G<:Braket.Gate, E<:OpenQASM3.Expression} = ctx(mod, argument, n_arg(argument), braket_gate) 
function (ctx::AbstractQASMContext)(mod::QuantumGateModifier{OpenQASM3.Inv, E}, argument, braket_gate::Vector{<:Instruction}) where {E<:OpenQASM3.Expression}
    to_map = reverse(braket_gate)
    braket_gate = map(to_map) do ix
        old_qc  = qubit_count(ix.operator)
        new_op  = inv(ix.operator)
        new_qc  = qubit_count(new_op)
        qc_diff = new_qc - old_qc
        new_ix  = Instruction(new_op, pad_qubits(new_op, qc_diff, ix.target))
    end
    return braket_gate
end
function (ctx::AbstractQASMContext)(mod::QuantumGateModifier{OpenQASM3.Pow, E}, argument, braket_gate::Vector{<:Instruction}) where {E<:OpenQASM3.Expression}
    braket_gate = map(braket_gate) do ix
        qc      = qubit_count(ix.operator)
        new_op  = ix.operator ^ argument
        new_ix  = Instruction(new_op, ix.target)
    end
    return braket_gate
end
(ctx::AbstractQASMContext)(mod::QuantumGateModifier{OpenQASM3.Inv, E}, argument, braket_gate::G) where {G<:Braket.Gate, E<:OpenQASM3.Expression} = inv(braket_gate)
(ctx::AbstractQASMContext)(mod::QuantumGateModifier{OpenQASM3.Pow, E}, argument, braket_gate::G) where {G<:Braket.Gate, E<:OpenQASM3.Expression} = braket_gate ^ argument
(ctx::AbstractQASMContext)(mod::QuantumGateModifier{OpenQASM3.NothingMod, OpenQASM3.NothingExpression}, argument, braket_gate) = braket_gate
function (ctx::AbstractQASMContext)(mods::Vector{M}, braket_gate) where {M<:OpenQASM3.QuantumGateModifier}
    for mod in Iterators.reverse(mods)
        braket_gate = ctx(mod, ctx(mod.argument), braket_gate)
    end
    return braket_gate
end
(ctx::AbstractQASMContext)(mods::Vector{QuantumGateModifier{OpenQASM3.NothingMod, OpenQASM3.NothingExpression}}, braket_gate) = braket_gate

_to_instructions(ix::Instruction{O}) where {O<:Braket.Operator} = [ix]
_to_instructions(ixs::Vector{<:Instruction}) = ixs
_lookup_gate(gate_name::String, gate_args::Vector{Float64}) = BuiltinGates[Symbol(gate_name)](gate_args...)
function _construct_raw_gate(ctx::AbstractQASMContext, gate_name::String, gate_args::Vector{Float64})
    sym_name = Symbol()
    try
        gd = lookup_def(GateDef, gate_name, ctx)
        return _to_instructions(gd(gate_args...))
    catch e
        return _lookup_gate(lowercase(gate_name), gate_args)
    end
end

function _splat_qubits(gate_qubits, qc::Int)
    max_length   = maximum(length, gate_qubits)
    final_qubits = Vector{Int}[Vector{Int}(undef, qc) for ix in 1:max_length]
    if !all(length(qubits) == max_length for qubits in gate_qubits)
        for r_ix in filter(r_ix->(length(gate_qubits[r_ix]) == 1), 1:length(gate_qubits))
            gate_qubits[r_ix] = [gate_qubits[r_ix][1] for _ in 1:max_length]
        end
    end
    for ii in 1:max_length, q in 1:qc
        final_qubits[ii][q] = gate_qubits[q][ii]
    end
    return final_qubits
end

function _attach_gate_to_qubits(ctx::AbstractQASMContext, braket_gate::Vector{<:Instruction}, qubits::Vector{I}) where {I<:OpenQASM3.AbstractIdentifier}
    qc          = qubit_count(braket_gate)
    gate_qubits = _get_qubits_from_mapping(qubits, ctx)
    if qc == 1
        for q in gate_qubits, ix in braket_gate
            push!(ctx, Instruction(ix.operator, q))
        end
    else
        length(gate_qubits) == qc || error("Gate def has qubit count $qc but $(length(gate_qubits)) were provided.")
        # handle splatting
        for qubits in _splat_qubits(gate_qubits, qc), ix in braket_gate
            push!(ctx, Instruction(ix.operator, qubits))
        end
    end
    return nothing
end

function _attach_gate_to_qubits(ctx::AbstractQASMContext, braket_gate::G, qubits::Vector{I}) where {I<:OpenQASM3.AbstractIdentifier, G<:Braket.Gate}
    qc          = qubit_count(braket_gate)
    gate_qubits = _get_qubits_from_mapping(qubits, ctx)
    if qc == 1
        for q in Iterators.flatten(gate_qubits)
            push!(ctx, Instruction{G}(braket_gate, q))
        end
    else
        length(gate_qubits) == qc || error("Gate of type $G has qubit count $qc but $(length(gate_qubits)) were provided.")
        # handle splatting
        for qubits in _splat_qubits(gate_qubits, qc) 
            push!(ctx, Instruction{G}(braket_gate, qubits))
        end
    end
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.QuantumGate{I, M, E}, gate_name::String) where {I, M, E}
    gate_args     = Float64[convert(Float64, ctx(a)) for a in node.arguments]
    gate          = _construct_raw_gate(ctx, gate_name, gate_args)
    modified_gate = ctx(node.modifiers, gate)
    _attach_gate_to_qubits(ctx, modified_gate, node.qubits)
    return
end

function (ctx::QASMSubroutineContext)(node::QuantumArgument, name::String, ix::Int)
    n_qubits = isnothing(node.size) ? 1 : ctx(node.size)
    ctx.definitions[name]   = QubitDef(n_qubits)
    ctx.qubit_mapping[name] = n_qubits == 1 ? [ix] : collect(ix:ix+n_qubits-1)
    for q_ix in 0:n_qubits - 1
        ctx.qubit_mapping[(name, q_ix)] = [ix]
        ix += 1
    end
    return ix
end

function (ctx::QASMSubroutineContext)(node::ClassicalArgument{T, AT}, name::String, ix::Int) where {T, AT}
    ctx.definitions[name] = ClassicalDef(undef, T, ClassicalVar()) 
    return ix
end

function _map_arguments_back(arg_passed::IE, arg_defined::OpenQASM3.ClassicalArgument{T, OpenQASM3.mutable}, fn_ctx) where {T, IE<:IndexExpression}
    outer_id   = id(arg_passed)
    inner_id   = id(arg_defined)
    ctx        = parent(fn_ctx)
    outer_arg  = ctx.definitions[outer_id]
    name, inds = _lookup_name(arg_passed, ctx)
    for (aix, ind) in enumerate(inds)
        ctx.definitions[name].value[ind.+1] = fn_ctx.definitions[inner_id].value[aix]
    end
    return
end
function _map_arguments_back(arg_passed::Identifier, arg_defined::OpenQASM3.ClassicalArgument{T, OpenQASM3.mutable}, fn_ctx) where {T}
    outer_id = id(arg_passed)
    inner_id = id(arg_defined)
    ctx      = parent(fn_ctx)
    ctx.definitions[outer_id] = fn_ctx.definitions[inner_id]
    return
end
_map_arguments_back(arg_passed, arg_defined::OpenQASM3.ClassicalArgument{T, OpenQASM3.readonly}, fn_ctx) where {T} = return
_map_arguments_back(arg_passed, arg_defined, fn_ctx) = return

function _map_arguments_forward(arg_passed, arg_defined::OpenQASM3.ClassicalArgument{T, AT}, qubit_alias::Dict, fn_ctx) where {T, AT}
    ctx       = parent(fn_ctx)
    arg_value = ctx(arg_passed)
    fn_ctx.definitions[id(arg_defined)] = ClassicalDef(arg_value, T, ClassicalVar())
    return
end

_get_alias_value(arg_val::String) = arg_val[1]
_get_alias_value(arg_val::Tuple{String,Vector{Int}}) = (arg_val[1], arg_val[2][1])
function _map_arguments_forward(arg_passed, arg_defined::OpenQASM3.QuantumArgument, qubit_alias::Dict, fn_ctx)
    ctx       = parent(fn_ctx)
    arg_name  = id(arg_defined)
    arg_val   = _lookup_name(arg_passed, ctx)
    qubit_alias[arg_name] = _get_alias_value(arg_val)
    qubit_alias[(arg_name, 0)] = _get_alias_value(arg_val)
    return
end

_to_literal(x::Int)  = OpenQASM3.IntegerLiteral(x)
_to_literal(x::Real) = OpenQASM3.FloatLiteral(Float64(x))
function _to_literal(x::BitVector)
    width = length(x)
    val = 0
    for (bi, b) in enumerate(x)
        val += b * 2^(bi-1)
    end
    return OpenQASM3.BitstringLiteral(val, width)
end
function (ctx::AbstractQASMContext)(node::OpenQASM3.FunctionCall)
    fn_name = id(node)
    if haskey(builtin_functions, fn_name)
        f = builtin_functions[fn_name]
        # wrap in OQ3 type
        args = map(_to_literal ∘ ctx, node.arguments)
        return f(args...)
    end
    fn_def  = lookup_def(SubroutineDef, fn_name, ctx; span=node.span)
    fn_ctx  = fn_def.ctx
    isnothing(fn_def) && error("Subroutine $fn_name not found!")
    qubit_alias = Dict()
    foreach(arg_pair->_map_arguments_forward(arg_pair[1], arg_pair[2], qubit_alias, fn_ctx), zip(node.arguments, fn_def.args))
    old_mapping           = deepcopy(fn_ctx.qubit_mapping)
    new_mapping           = Dict(kq=>ctx.qubit_mapping[qubit_alias[kq]] for (kq, vq) in fn_ctx.qubit_mapping)
    fn_ctx.qubit_mapping  = new_mapping
    rv, out               = fn_def.body_gen()
    fn_ctx.qubit_mapping  = old_mapping
    # map arguments back to outer context
    foreach(arg_pair->_map_arguments_back(arg_pair[1], arg_pair[2], fn_ctx), zip(node.arguments, fn_def.args))
    return rv
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.SubroutineDefinition{S}, subroutine_name::String) where {S}
    subroutine_body = node.body
    subroutine_args = node.arguments
    block_ctx       = _new_subroutine_scope(ctx)
    ix = 0
    for arg in subroutine_args 
        ix = block_ctx(arg, ix)
    end
    function subroutine_body_builder()
        return_value    = nothing
        for statement in subroutine_body
            block_ctx(statement)
            if statement isa ReturnStatement
                return_value = block_ctx(statement.expression)
                break
            end
        end
        return return_value, output(block_ctx)
    end
    ctx.definitions[subroutine_name] = SubroutineDef(node.arguments, subroutine_body_builder, node.return_type, block_ctx)
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.ForInLoop)
    _check_annotations(node)
    loopvar_name  = id(node)
    loop_set      = ctx(node.set_declaration)
    block_ctx     = _new_inner_scope(ctx)
    for i in loop_set
        scalar_matches_type(node.type, i, "Loop variable type error: $i is not a ")
        block_ctx.definitions[loopvar_name] = ClassicalDef(i, node.type, ClassicalVar())
        for subnode in node.block
            maybe_break_or_continue = block_ctx(subnode)
            if maybe_break_or_continue isa OpenQASM3.BreakStatement
                return nothing
            elseif maybe_break_or_continue isa OpenQASM3.ContinueStatement
                break
            end
        end
    end
    # update variables in the parent defs
    for k in intersect(keys(ctx.definitions), keys(block_ctx.definitions))
        ctx.definitions[k] = block_ctx.definitions[k]
    end
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.WhileLoop)
    _check_annotations(node)
    block_ctx = _new_inner_scope(ctx)
    while block_ctx(node.while_condition)
        for subnode in node.block
            maybe_break_or_continue = block_ctx(subnode)
            if maybe_break_or_continue isa OpenQASM3.BreakStatement
                return nothing
            elseif maybe_break_or_continue isa OpenQASM3.ContinueStatement
                break
            end
        end
    end
    # update variables in the parent defs
    for k in intersect(keys(ctx.definitions), keys(block_ctx.definitions))
        ctx.definitions[k] = block_ctx.definitions[k]
    end
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.BranchingStatement)
    _check_annotations(node)
    block_ctx = _new_inner_scope(ctx)
    if ctx(node.condition) > 0
        for subnode in node.if_block
            res = block_ctx(subnode)
            isbreakorcontinue(res) && return res
        end
    elseif !isnothing(node.else_block)
        for subnode in node.else_block
            res = block_ctx(subnode)
            isbreakorcontinue(res) && return res
        end
    end
    return nothing
end

function _check_undefined(node, name, ctx; local_only=true)
    is_defined(name, ctx; local_only) && error("Identifier $name already in use at $(node.span).")
    return nothing
end
_check_undefined(node::ClassicalAssignment, name, ctx; local_only=true) = return
_check_undefined(node::ClassicalDeclaration, name, ctx; local_only=true) = return
_check_undefined(node::IndexedIdentifier, name, ctx; local_only=true) = return
_check_undefined(node::QuantumGate, name, ctx; local_only=true) = return

isbreakorcontinue(node::OpenQASM3.BreakStatement)    = true
isbreakorcontinue(node::OpenQASM3.ContinueStatement) = true
isbreakorcontinue(node::Nothing) = false
isbreakorcontinue(node)          = false

_is_hardware_qubit(qname::AbstractString) = !isnothing(match(r"\$[0-9]+", qname))
_is_hardware_qubit(qname::OpenQASM3.Identifier) = _is_hardware_qubit(id(qname))
_is_hardware_qubit(qname) = false
resolve_qubit(q::Tuple{String, Int}, q_name::String, def::QubitDef, ctx::AbstractQASMContext)   = [q]
resolve_qubit(q::String, q_name::String, def::QubitDef, ctx::AbstractQASMContext)               = [(q, ix) for ix in 0:def.size-1]
resolve_qubit(q::OpenQASM3.Identifier, q_name::String, def::QubitDef, ctx::AbstractQASMContext) = resolve_qubit(q_name, q_name, def, ctx) 
resolve_qubit(q::Tuple{String, T}, q_name::String, def::QubitDef, ctx::AbstractQASMContext) where {T} = [(q_name, ix) for ix in Iterators.flatten(q[2])]
_get_indices(ix::Int)::Vector{Int} = [ix]
_get_indices(ix::Vector{Int})::Vector{Int} = ix
_get_indices(ix::Vector{Vector{Int}})::Vector{Int} = ix[1]
function resolve_qubit(node::OpenQASM3.IndexedIdentifier{IT}, q_name::String, def::QubitDef, ctx::AbstractQASMContext) where {IT<:OpenQASM3.IndexElement}
    ixs = _get_indices(ctx(node.indices))
    return resolve_qubit((q_name, ixs), q_name, def, ctx) 
end
function resolve_qubit(q, q_name::String, ctx::AbstractQASMContext, span=())
    _is_hardware_qubit(q_name) && return [id(q_name)]
    def = lookup_def(QubitDef, q_name, ctx; span=span)
    return resolve_qubit(q, q_name, def, ctx) 
end
resolve_qubit(q::String, ctx::AbstractQASMContext, span=()) = _is_hardware_qubit(q) ? [id(q)] : resolve_qubit(q, id(q), ctx, span)
resolve_qubit(q, ctx::AbstractQASMContext, span=()) = resolve_qubit(q, id(q), ctx, span)

function _get_qubits_from_mapping(q::Vector{I}, ctx::AbstractQASMContext) where {I}
    qubits = Vector{Vector{Int}}(undef, length(q))
    for (i, qubit) in enumerate(q)
        qubits[i] = Int[]
        resolved_q = resolve_qubit(qubit, ctx)
        for rq in resolved_q
            if _is_hardware_qubit(rq)
                push!(qubits[i], tryparse(Int, replace(rq, "\$"=>"")))
                continue
            else
                append!(qubits[i], Iterators.flatten(ctx.qubit_mapping[rq]))
            end
        end
    end
    return qubits
end
_get_qubits_from_mapping(q::Tuple{String, Int}, ctx::AbstractQASMContext) = ctx.qubit_mapping[q]
# does nothing for now
function (ctx::AbstractQASMContext)(node::Union{OpenQASM3.QuantumMeasurement, OpenQASM3.QuantumMeasurementStatement})
    _check_annotations(node)
    return nothing
end

function (ctx::QASMSubroutineContext)(node::OpenQASM3.ReturnStatement)
    _check_annotations(node)
    return ctx(node.expression)
end

_process_init_expression(::Nothing, type::T, name::String, ctx::AbstractQASMContext) where {T} = nothing 
_process_init_expression(::OpenQASM3.NothingExpression, type::T, name::String, ctx::AbstractQASMContext) where {T} = nothing
_process_init_expression(expr::OpenQASM3.Expression, type::T, name::String, ctx::AbstractQASMContext) where {T} = ctx(expr)
for T in (:(OpenQASM3.FloatType), :(OpenQASM3.ComplexType), :(OpenQASM3.IntType), :(OpenQASM3.UintType), :(OpenQASM3.BoolType), :(OpenQASM3.BitType))
    @eval begin
        value_type(::Type{$T}, ::OpenQASM3.NothingExpression) = Nothing
    end
end
value_type(::Type{OpenQASM3.FloatType}, exp::OpenQASM3.Expression) = Float64
value_type(::Type{OpenQASM3.ComplexType}, exp::OpenQASM3.Expression) = ComplexF64
value_type(::Type{OpenQASM3.IntType}, exp::OpenQASM3.Expression) = Int
value_type(::Type{OpenQASM3.UintType}, exp::OpenQASM3.Expression) = Int
value_type(::Type{OpenQASM3.BoolType}, exp::OpenQASM3.Expression) = Bool
value_type(::Type{OpenQASM3.BitType}, exp::OpenQASM3.Expression) = BitVector
value_type(::Type{OpenQASM3.ArrayType{T}}, exp::OpenQASM3.Expression) where {T} = Vector{value_type(T, exp)}
function (ctx::AbstractQASMContext)(node::OpenQASM3.ClassicalDeclaration{T, IE}, name::String) where {T<:ClassicalType, IE}
    val_T = value_type(T, node.init_expression)
    val   = _process_init_expression(node.init_expression, node.type, name, ctx)
    ctx.definitions[name] = ClassicalDef{val_T, T, ClassicalVar}(val, node.type, ClassicalVar())
    return nothing
end

function (ctx::AbstractQASMContext)(node::OpenQASM3.ConstantDeclaration{T}, name::String) where {T}
    val  = _process_init_expression(node.init_expression, T, name, ctx)
    ctx.definitions[name] = ClassicalDef(val, T, ClassicalConst())
    return nothing
end

_check_lval(lv::OpenQASM3.Identifier)        = return
_check_lval(lv::OpenQASM3.IndexedIdentifier) = return
_check_lval(lv::T) where {T} = error("Assignment not implemented for $T.")
_check_node_op(op::OpenQASM3.AssignmentOperator) = return
_check_node_op(op) = error("Unknown op $op in assigment.")

_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(=)}, ctx::AbstractQASMContext) = r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(+=)}, ctx::AbstractQASMContext) = d_value + r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(-=)}, ctx::AbstractQASMContext) = d_value - r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(*=)}, ctx::AbstractQASMContext) = d_value * r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(/=)}, ctx::AbstractQASMContext) = d_value / r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(|=)}, ctx::AbstractQASMContext) = d_value | r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(&=)}, ctx::AbstractQASMContext) = d_value & r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{:(^=)}, ctx::AbstractQASMContext) = d_value ⊻ r_val
_assign_lvalue(lvalue::Identifier, r_val::Number, d_value::Number, ::Val{O}, ctx::AbstractQASMContext) where {O} = error("Operation $op not defined.")
_assign_lvalue(lvalue::Identifier, r_val::Vector, d_value::Number, op, ctx::AbstractQASMContext) = _assign_lvalue(lvalue, r_val[1], d_value, op, ctx)
function _assign_lvalue(lvalue::IndexedIdentifier, r_val::AbstractVector, d_value::AbstractVector, op, ctx::AbstractQASMContext)
    len = length(lookup_def(QASMDefinition, id(lvalue), ctx))
    l_inds = collect(Iterators.flatten(mapreduce(ix->ctx(ix, len), vcat, lvalue.indices)))
    for (ii, l_ind) in enumerate(l_inds)
        d_value[l_ind+1] = _assign_lvalue(lvalue.name, r_val[ii], d_value[l_ind+1], op, ctx)
    end
    return d_value
end
_assign_lvalue(lvalue::IndexedIdentifier, r_val::Number, d_value::AbstractVector, op, ctx::AbstractQASMContext) = _assign_lvalue(lvalue, fill(r_val, length(d_value)), d_value, op, ctx)
_assign_lvalue(lvalue::IndexedIdentifier, r_val::Number, d_value::Number, op, ctx::AbstractQASMContext) = _assign_lvalue(lvalue, [r_val], [d_value], op, ctx)
_assign_lvalue(lvalue, r_val, d_value::Nothing, op, ctx::AbstractQASMContext) = r_val

function (ctx::AbstractQASMContext)(node::OpenQASM3.ClassicalAssignment{O}, name::String) where {O}
    _check_lval(node.lvalue) 
    def  = lookup_def(ClassicalDef, name, ctx; span=node.span)
    def.kind isa ClassicalVar || error("Variable `$name` cannot be assigned to.")
    val   = _assign_lvalue(node.lvalue, ctx(node.rvalue), def.value, O(), ctx)
    ctx.definitions[name] = ClassicalDef(val, def.valtype, def.kind)
    return
end

function _add_to_qubit_mapping(ctx::QASMGlobalContext, num_qubits::Int, qname::String)
    ctx.qubit_mapping[qname] = collect(ctx.n_qubits:(ctx.n_qubits + num_qubits - 1))
    for q_ix in 0:num_qubits - 1
        ctx.qubit_mapping[(qname, q_ix)] = [ctx.n_qubits + q_ix]
    end
    ctx.n_qubits += num_qubits
    return
end
_add_to_qubit_mapping(ctx::QASMGlobalContext, ::Nothing, qname::String) = _add_to_qubit_mapping(ctx, 1, qname)
function (ctx::QASMGlobalContext)(node::OpenQASM3.QubitDeclaration, qname::String)
    num_qubits               = value(node.size)
    ctx.definitions[qname]   = QubitDef(num_qubits)
    _add_to_qubit_mapping(ctx, num_qubits, qname)
    return nothing
end

function (ctx::QASMGlobalContext)(node::OpenQASM3.Program)
    for s in node.statements
        res = ctx(s)
        isbreakorcontinue(res) && ctx(res)
    end
end


function Base.collect(ctx::QASMGlobalContext)
    c = Circuit()
    wo = output(ctx)
    foreach(ix->Braket.add_instruction!(c, ix), wo.ixs)
    for rt in wo.results
        obs = Braket.extract_observable(rt)
        if !isnothing(obs) && c.observables_simultaneously_measureable && !(rt isa AdjointGradient)
            Braket.add_to_qubit_observable_mapping!(c, obs, rt.targets)
        end
        Braket.add_to_qubit_observable_set!(c, rt)
        push!(c.result_types, rt)
    end
    return c
end
function interpret(program::OpenQASM3.Program, extern_lookup::Dict{String, <:Any}=Dict{String,Float64}())
    global_ctx = QASMGlobalContext{Braket.Operator}(extern_lookup)
    # walk the nodes recursively
    global_ctx(program)
    return collect(global_ctx)
end
interpret(program::OpenQASM3.Program, ::Nothing) = interpret(program, Dict{String,Float64}())
