module Quasar

using Automa, AbstractTrees, DataStructures, Dates
using DataStructures: Stack

export parse_qasm, QasmProgramVisitor

struct QasmParseError <: Exception
    message::String
    parse_stack::Stack
    position::Int
    qasm::String
end
function Base.showerror(io::IO, err::QasmParseError)
    print(io, "QasmParseError: ")
    print(io, err.message * "\n")
    max_codeunits = min(length(err.qasm), err.position+100)
    print(io, "Qasm location: ", err.qasm[err.position:max_codeunits])
end

include("builtin_functions.jl")
#const unicode = Automa.RE(read(joinpath(@__DIR__, "unicode_re.dat"), String))
const unicode = re"α|β|γ|δ|ϵ|ε|ϕ|φ|ζ|η|θ|ϑ|ι|κ|λ|μ|ν|ξ|ρ|σ|ς|υ|χ|ψ|ω"
const first_letter   = re"[A-Za-z_]" | unicode
const general_letter = first_letter | re"[0-9]" 

const prefloat = re"[-+]?([0-9]+\.[0-9]*|[0-9]*\.[0-9]+)"
const integer = re"[-+]?[0-9]+"
const float   = prefloat | ((prefloat | re"[-+]?[0-9]+") * re"[eE][-+]?[0-9]+")

const qasm_tokens = [
        :identifier   => first_letter * rep(general_letter),
        :irrational   => re"π|pi|τ|tau|ℯ|ℇ|euler",
        :comma        => re",",
        :colon        => re":",
        :semicolon    => re";",
        :question     => re"\?",
        :equal        => re"=",
        :lparen       => re"\(",
        :rparen       => re"\)",
        :lbracket     => re"\[",
        :rbracket     => re"]",
        :lbrace       => re"{",
        :rbrace       => re"}",
        :annot        => re"@[*]",
        :at           => re"@",
        :version      => re"OPENQASM",
        :input        => re"input",
        :output       => re"output",
        :pragma       => re"#pragma",
        :qubit        => re"qubit",
        :hw_qubit     => re"$[0-9]+",
        :gate_def     => re"gate",
        :function_def => re"def",
        :if_block     => re"if",
        :else_block   => re"else",
        :switch_block => re"switch",
        :while_block  => re"while",
        :in_token     => re"in",
        :for_block    => re"for",
        :return_token => re"return",
        :control_mod  => re"ctrl",
        :negctrl_mod  => re"negctrl",
        :inverse_mod  => re"inv",
        :power_mod    => re"pow",
        :measure      => re"measure",
        :arrow_token  => re"->",
        :reset_token  => re"reset",
        :delay_token  => re"delay",
        :barrier_token => re"barrier",
        :void         => re"void",
        :const_token  => re"const",
        :assignment   => re"=|-=|\+=|\*=|/=|^=|&=|\|=|<<=|>>=",
        :operator     => re"-|\+|\++|\*|\*\*|/|%|&|&&|\||\|\||^|!|!=|~|>|<|<<|>>|>=|<=|=>|==",
        :boolean      => re"true|false",
        :bitstring    => re"\"([01] _?)* [01]\"",
        :all_token    => re"all",
        :break_token  => re"break",
        :mutable      => re"mutable",
        :readonly     => re"readonly",
        :builtin_gate => re"gphase|U",
        :alias        => re"let",
        :box          => re"box",
        :end_token    => re"end",
        :dim_token    => re"#dim[ ]?=[ ]?[0-7]",
        :im_token     => re"im",
        :case         => re"case",
        :default      => re"default",
        :keyword      => re"creg|qreg",
        :oct          => re"0o[0-7]+",
        :bin          => re"(0b|0B)[0-1]+",
        :hex          => re"0x[0-9A-Fa-f]+",
        :dot          => re"\.",
        :integer_token => integer,
        :float_token   => float,
        :include_token   => re"include",
        :continue_token  => re"continue",
        :octal_integer   => re"0o([0-7]_?)* [0-7]",
        :hex_integer     => re"(0x|0X) ([0-9a-fA-F] _?)* [0-9a-fA-F]",
        :hardware_qubit  => re"$ [0-9]+",
        :line_comment    => re"//",
        :block_comment   => re"/\* .*? \*/",
        :char            => '\'' * (re"[ -&(-~]" | ('\\' * re"[ -~]")) * '\'',
        :string_token    => '"' * rep(re"[ !#-~]" | re"\\\\\"") * '"' | '\'' * rep(re"[ -&(-~]" | ('\\' * re"[ -~]")) * '\'',
        :newline         => re"\r?\n",
        :spaces          => re"[\t ]+",
        :classical_type    => re"bool|uint|int|float|angle|complex|array|bit|stretch|duration",
        :durationof_token  => re"durationof", # this MUST be lower than classical_type to preempt duration
        :duration_literal  => (float | integer) * re"dt|ns|us|ms|s|\xce\xbc\x73", # transcode'd μs
        :forbidden_keyword => re"cal|defcal|extern",
]

const dt_type = Ref{DataType}()

function __init__()
    dt_type[] = Nanosecond
end
function parse_pragma end
function visit_pragma end

@eval @enum Token error $(first.(qasm_tokens)...)
make_tokenizer((error,
    [Token(i) => j for (i,j) in enumerate(last.(qasm_tokens))]
)) |> eval

struct QasmExpression
    head::Symbol
    args::Vector{Any}
    QasmExpression(head::Symbol, args::Vector) = new(head, args)
end
QasmExpression(head) = QasmExpression(head, [])
QasmExpression(head, @nospecialize(args...)) = QasmExpression(head, collect(args))
QasmExpression(head, arg) = QasmExpression(head, [arg])

Base.show(io::IO, qasm_expr::QasmExpression) = print_tree(io, qasm_expr, maxdepth=10)
Base.iterate(qasm_expr::QasmExpression) = (qasm_expr, nothing)
Base.iterate(qasm_expr::QasmExpression, ::Nothing) = nothing
Base.length(qasm_expr::QasmExpression) = 1
Base.push!(qasm_expr::QasmExpression, arg) = push!(qasm_expr.args, arg)
Base.append!(qasm_expr::QasmExpression, arg::QasmExpression) = push!(qasm_expr.args, arg)
Base.append!(qasm_expr::QasmExpression, args::Vector) = append!(qasm_expr.args, args)
Base.pop!(qasm_expr::QasmExpression) = pop!(qasm_expr.args)
Base.copy(qasm_expr::QasmExpression) = QasmExpression(qasm_expr.head, deepcopy(qasm_expr.args))

head(qasm_expr::QasmExpression) = qasm_expr.head

AbstractTrees.children(qasm_expr::QasmExpression) = qasm_expr.args
AbstractTrees.printnode(io::IO, qasm_expr::QasmExpression) = print(io, "QasmExpression :$(qasm_expr.head)")

function Base.:(==)(qasm_a::QasmExpression, qasm_b::QasmExpression)
    a_children = children(qasm_a)
    b_children = children(qasm_b)
    length(a_children) != length(b_children) && return false
    return a_children == b_children
end

parse_hw_qubit(token, qasm)   = QasmExpression(:hw_qubit, qasm[token[1]:token[1]+token[2]-1])
parse_identifier(token, qasm) = QasmExpression(:identifier, String(codeunits(qasm)[token[1]:token[1]+token[2]-1]))
function extract_scope(tokens, stack, start, qasm)
    # a "scope" begins with an { and ends with an }
    # but we may have nested scope!
    opener = popfirst!(tokens)
    opener[end] == lbrace || throw(QasmParseError("scope does not open with {", stack, start, qasm))
    # need to match openers to closers to exit the scope
    openers_met  = 1
    closers_met  = 0
    scope_tokens = Tuple{Int64, Int32, Token}[]
    while closers_met < openers_met && !isempty(tokens)
        next_token      = popfirst!(tokens)
        next_token[end] == lbrace && (openers_met += 1)
        next_token[end] == rbrace && (closers_met += 1)
        push!(scope_tokens, next_token)
    end
    pop!(scope_tokens) # closing }
    return scope_tokens
end

function parse_scope(tokens, stack, start, qasm)
    scope_tokens = extract_scope(tokens, stack, start, qasm)
    return parse_qasm(scope_tokens, qasm, QasmExpression(:scope))
end

function parse_block_body(expr, tokens, stack, start, qasm)
    is_scope = tokens[1][end] == lbrace
    if is_scope
        body       = parse_scope(tokens, stack, start, qasm)
        body_exprs = convert(Vector{QasmExpression}, collect(Iterators.reverse(body)))::Vector{QasmExpression}
        foreach(body_expr->push!(body_exprs[1], body_expr), body_exprs[2:end])
        push!(expr, body_exprs[1])
    else # one line
        eol = findfirst(triplet->triplet[end] == semicolon, tokens)
        body_tokens = splice!(tokens, 1:eol)
        body = parse_expression(body_tokens, stack, start, qasm)
        push!(expr, body)
    end
end

function parse_switch_block(tokens, stack, start, qasm)
    expr = QasmExpression(:switch)
    cond_tokens = extract_parensed(tokens, stack, start, qasm)
    push!(cond_tokens, (-1, Int32(-1), semicolon))
    push!(expr, parse_expression(cond_tokens, stack, start, qasm))
    interior_tokens = extract_scope(tokens, stack, start, qasm)
    met_default = false
    while !isempty(interior_tokens)
        next_token = popfirst!(interior_tokens)
        if next_token[end] == case
            !met_default || throw(QasmParseError("case statement cannot occur after default in switch block.", stack, start, qasm))
            case_expr = QasmExpression(:case)
            brace_loc = findfirst(triplet->triplet[end] == lbrace, interior_tokens)
            isnothing(brace_loc) && throw(QasmParseError("case statement missing opening {", stack, start, qasm))
            val_tokens = splice!(interior_tokens, 1:brace_loc-1)
            push!(val_tokens, (-1, Int32(-1), semicolon))
            push!(case_expr, parse_list_expression(val_tokens, stack, start, qasm))
            parse_block_body(case_expr, interior_tokens, stack, start, qasm)
            push!(expr, case_expr)
        elseif next_token[end] == default
            !met_default || throw(QasmParseError("only one default statement allowed in switch block.", stack, start, qasm))
            default_expr = QasmExpression(:default)
            parse_block_body(default_expr, interior_tokens, stack, start, qasm)
            push!(expr, default_expr)
            met_default = true
        elseif next_token[end] == newline
            continue
        else
            throw(QasmParseError("invalid switch-case statement.", stack, start, qasm))
        end
    end
    return expr
end

function parse_if_block(tokens, stack, start, qasm)
    condition_tokens = extract_parensed(tokens, stack, start, qasm)
    push!(condition_tokens, (-1, Int32(-1), semicolon))
    condition_value = parse_expression(condition_tokens, stack, start, qasm)
    if_expr = QasmExpression(:if, condition_value)
    # handle condition
    parse_block_body(if_expr, tokens, stack, start, qasm)
    has_else = tokens[1][end] == else_block 
    if has_else
        popfirst!(tokens) # else
        else_expr = QasmExpression(:else)
        parse_block_body(else_expr, tokens, stack, start, qasm)
        push!(if_expr, else_expr)
    end
    return if_expr
end
function parse_while_loop(tokens, stack, start, qasm)
    condition_tokens = extract_parensed(tokens, stack, start, qasm)
    push!(condition_tokens, (-1, Int32(-1), semicolon))
    condition_value = parse_expression(condition_tokens, stack, start, qasm)
    while_expr = QasmExpression(:while, condition_value)
    # handle condition
    parse_block_body(while_expr, tokens, stack, start, qasm)
    return while_expr
end
function parse_for_loop(tokens, loop_var_type, loop_var_name, loop_values, stack, start, qasm)
    for_expr = QasmExpression(:for, loop_var_type, loop_var_name, loop_values)
    parse_block_body(for_expr, tokens, stack, start, qasm)
    return for_expr
end

function parse_arguments_list(tokens, stack, start, qasm)
    arguments = QasmExpression(:arguments)
    first(tokens)[end] != lparen && return arguments
    interior = extract_parensed(tokens, stack, start, qasm)
    push!(interior, (-1, Int32(-1), semicolon))
    push!(arguments, parse_list_expression(interior, stack, start, qasm))
    return arguments
end

function parse_function_def(tokens, stack, start, qasm)
    function_name    = popfirst!(tokens)
    function_name[end] == identifier || throw(QasmParseError("function definition must have a valid identifier as a name", stack, start, qasm))
    function_name_id = parse_identifier(function_name, qasm) 
    arguments        = parse_arguments_list(tokens, stack, start, qasm)
    has_return_type  = tokens[1][end] == arrow_token
    if has_return_type
        arrow = popfirst!(tokens)
        tokens[1][end] == classical_type || throw(QasmParseError("function return type must be a classical type", stack, start, qasm))
        return_type = parse_classical_type(tokens, stack, start, qasm)
    else
        return_type = QasmExpression(:void)
    end
    expr = QasmExpression(:function_definition, function_name_id, arguments, return_type)
    parse_block_body(expr, tokens, stack, start, qasm)
    return expr
end
function parse_gate_def(tokens, stack, start, qasm)
    gate_name = popfirst!(tokens)
    gate_name[end] == identifier || throw(QasmParseError("gate definition must have a valid identifier as a name", stack, start, qasm))
    gate_name_id = parse_identifier(gate_name, qasm)

    gate_args    = parse_arguments_list(tokens, stack, start, qasm)
    qubit_tokens = splice!(tokens, 1:findfirst(triplet->triplet[end]==lbrace, tokens)-1)
    push!(qubit_tokens, (-1, Int32(-1), semicolon))
    target_expr = QasmExpression(:qubit_targets, parse_list_expression(qubit_tokens, stack, start, qasm))
    expr = QasmExpression(:gate_definition, gate_name_id, gate_args, target_expr)
    parse_block_body(expr, tokens, stack, start, qasm)
    return expr
end

struct SizedBitVector <: AbstractArray{Bool, 1}
    size::QasmExpression
    SizedBitVector(size::QasmExpression) = new(size)
    SizedBitVector(sbv::SizedBitVector) = new(sbv.size)
end
Base.length(s::SizedBitVector) = s.size
Base.size(s::SizedBitVector) = (s.size,)
Base.show(io::IO, s::SizedBitVector) = print(io, "SizedBitVector{$(s.size.args[end])}")
struct SizedInt <: Integer
    size::QasmExpression
    SizedInt(size::QasmExpression) = new(size)
    SizedInt(sint::SizedInt) = new(sint.size)
end
Base.show(io::IO, s::SizedInt) = print(io, "SizedInt{$(s.size.args[end])}")
struct SizedUInt <: Unsigned 
    size::QasmExpression
    SizedUInt(size::QasmExpression) = new(size)
    SizedUInt(suint::SizedUInt) = new(suint.size)
end
Base.show(io::IO, s::SizedUInt) = print(io, "SizedUInt{$(s.size.args[end])}")
struct SizedFloat <: AbstractFloat
    size::QasmExpression
    SizedFloat(size::QasmExpression) = new(size)
    SizedFloat(sfloat::SizedFloat) = new(sfloat.size)
end
Base.show(io::IO, s::SizedFloat) = print(io, "SizedFloat{$(s.size.args[end])}")
struct SizedAngle <: AbstractFloat
    size::QasmExpression
    SizedAngle(size::QasmExpression) = new(size)
    SizedAngle(sangle::SizedAngle) = new(sangle.size)
end
Base.show(io::IO, s::SizedAngle) = print(io, "SizedAngle{$(s.size.args[end])}")
struct SizedComplex <: Number
    size::QasmExpression
    SizedComplex(size::QasmExpression) = new(size)
    SizedComplex(scomplex::SizedComplex) = new(scomplex.size)
end
Base.show(io::IO, s::SizedComplex) = print(io, "SizedComplex{$(s.size.args[end])}")

struct SizedArray{T,N} <: AbstractArray{T, N} 
    type::T
    size::NTuple{N, Int}
end
function SizedArray(eltype::QasmExpression, size::QasmExpression)
    arr_size = if head(size) == :n_dims
        ntuple(i->0, size.args[1].args[1])
    else
        ntuple(i->size.args[i], length(size.args))
    end
    return SizedArray(eltype.args[1], arr_size)
end
Base.show(io::IO, s::SizedArray{T, N}) where {T, N} = print(io, "SizedArray{$(sprint(show, s.type)), $N}")
Base.size(a::SizedArray{T, N}, dim::Int=0) where {T, N} = a.size[dim+1]

const SizedNumber = Union{SizedComplex, SizedAngle, SizedFloat, SizedInt, SizedUInt}

function parse_classical_type(tokens, stack, start, qasm)
    is_sized  = length(tokens) > 1 && tokens[2][end] == lbracket
    type_name = popfirst!(tokens)
    type_name[end] == classical_type || throw(QasmParseError("classical variable must have a classical type", stack, start, qasm))
    var_type = qasm[type_name[1]:type_name[1]+type_name[2]-1]
    if var_type == "complex"
        complex_tokens = extract_braced_block(tokens, stack, start, qasm)
        eltype = parse_classical_type(complex_tokens, stack, start, qasm)
        size   = eltype.args[1].size
    elseif var_type == "array"
        array_tokens = extract_braced_block(tokens, stack, start, qasm)
        eltype = parse_classical_type(array_tokens, stack, start, qasm)
        first(array_tokens)[end] == comma && popfirst!(array_tokens)
        size   = parse_expression(array_tokens, stack, start, qasm)
        return QasmExpression(:classical_type, SizedArray(eltype, size))
    elseif var_type == "duration"
        @warn "duration expression encountered -- currently `duration` is a no-op"
        # TODO: add proper parsing of duration expressions, including
        # support for units and algebraic durations like 2*a.
        return QasmExpression(:classical_type, :duration)
    elseif var_type == "stretch"
        @warn "stretch expression encountered -- currently `stretch` is a no-op"
        # TODO: add proper parsing of stretch expressions
        return QasmExpression(:classical_type, :stretch)
    else
        !any(triplet->triplet[end] == semicolon, tokens) && push!(tokens, (-1, Int32(-1), semicolon))
        size = is_sized ? parse_expression(tokens, stack, start, qasm) : QasmExpression(:integer_literal, -1)

    end
    if var_type == "bit"
        return QasmExpression(:classical_type, SizedBitVector(size))
    elseif var_type == "int"
        return QasmExpression(:classical_type, SizedInt(size))
    elseif var_type == "uint"
        return QasmExpression(:classical_type, SizedUInt(size))
    elseif var_type == "float"
        return QasmExpression(:classical_type, SizedFloat(size))
    elseif var_type == "complex"
        return QasmExpression(:classical_type, SizedComplex(size))
    elseif var_type == "bool"
        return QasmExpression(:classical_type, Bool)
    end
    throw(QasmParseError("could not parse classical type", stack, start, qasm))
end

function parse_classical_var(tokens, stack, start, qasm)
    # detect if we have a declared size
    is_declaration = false
    if tokens[end][end] == semicolon
        is_declaration = true
        pop!(tokens)
    end
    name = pop!(tokens)
    type = parse_classical_type(tokens, stack, start, qasm)
    name[end] == identifier || throw(QasmParseError("classical variable must have a valid name", stack, start, qasm))
    return type, parse_identifier(name, qasm)
end

const binary_assignment_ops = Dict{String, Symbol}(
                                                   "="=>Symbol("="),
                                                   "-="=>Symbol("-"),
                                                   "+="=>Symbol("+"),
                                                   "*="=>Symbol("*"),
                                                   "/="=>Symbol("/"),
                                                   "^="=>Symbol("^"),
                                                   "&="=>Symbol("&"),
                                                   "|="=>Symbol("|"),
                                                   "<<="=>Symbol("<<"),
                                                   ">>="=>Symbol(">>"),
                                                  )
function parse_assignment_op(op_token, qasm)
    op_string = parse_identifier(op_token, qasm)
    return binary_assignment_ops[op_string.args[1]::String]
end

parse_string_literal(token, qasm)  = QasmExpression(:string_literal, String(qasm[token[1]:token[1]+token[2]-1]))
parse_integer_literal(token, qasm) = QasmExpression(:integer_literal, tryparse(Int, qasm[token[1]:token[1]+token[2]-1]))
parse_hex_literal(token, qasm)     = QasmExpression(:integer_literal, tryparse(UInt, qasm[token[1]:token[1]+token[2]-1]))
parse_oct_literal(token, qasm)     = QasmExpression(:integer_literal, tryparse(Int, qasm[token[1]:token[1]+token[2]-1]))
parse_bin_literal(token, qasm)     = QasmExpression(:integer_literal, tryparse(Int, qasm[token[1]:token[1]+token[2]-1]))
parse_float_literal(token, qasm)   = QasmExpression(:float_literal, tryparse(Float64, qasm[token[1]:token[1]+token[2]-1]))
parse_boolean_literal(token, qasm)   = QasmExpression(:boolean_literal, tryparse(Bool, qasm[token[1]:token[1]+token[2]-1]))
function parse_duration_literal(token, qasm)
    str = String(codeunits(qasm)[token[1]:token[1]+token[2]-1])
    duration = if endswith(str, "ns")
            Nanosecond(tryparse(Int, chop(str, tail=2)))
        elseif endswith(str, "ms")
            Millisecond(tryparse(Int, chop(str, tail=2)))
        elseif endswith(str, "us") || endswith(str, "μs")
            Microsecond(tryparse(Int, chop(str, tail=2)))
        elseif endswith(str, "s")
            Second(tryparse(Int, chop(str, tail=1)))
        elseif endswith(str, "dt")
            dt_type[](tryparse(Int, chop(str, tail=2)))
        end
    QasmExpression(:duration_literal, duration)
end
function parse_irrational_literal(token, qasm)
    raw_string = String(codeunits(qasm)[token[1]:token[1]+token[2]-1])
    raw_string == "pi" && return QasmExpression(:irrational_literal, π)
    raw_string == "euler" && return QasmExpression(:irrational_literal, ℯ)
    raw_string == "tau" && return QasmExpression(:irrational_literal, 2*π)
    raw_string == "π" && return QasmExpression(:irrational_literal, π)
    raw_string == "τ" && return QasmExpression(:irrational_literal, 2*π)
    raw_string ∈ ("ℯ", "ℇ") && return QasmExpression(:irrational_literal, ℯ)
end
function parse_set_expression(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    interior = extract_scope(tokens, stack, start, qasm)
    set_elements = QasmExpression(:array_literal)
    push!(interior, (-1, Int32(-1), semicolon))
    while !isempty(interior)
        push!(set_elements, parse_expression(interior, stack, start, qasm))
        next_token = first(interior)
        next_token[end] == comma && popfirst!(interior)
        next_token[end] == semicolon && break 
    end
    return set_elements
end

function extract_braced_block(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    bracket_opening = findfirst(triplet->triplet[end] == lbracket, tokens)
    bracket_closing = findfirst(triplet->triplet[end] == rbracket, tokens)
    isnothing(bracket_opening) && throw(QasmParseError("missing opening [ ", stack, start, qasm))
    opener = popat!(tokens, bracket_opening)
    openers_met  = 1
    closers_met  = 0
    braced_tokens = Tuple{Int64, Int32, Token}[]
    while closers_met < openers_met && !isempty(tokens)
        next_token      = popfirst!(tokens)
        next_token[end] == lbracket && (openers_met += 1)
        next_token[end] == rbracket && (closers_met += 1)
        push!(braced_tokens, next_token)
    end
    pop!(braced_tokens) # closing ]
    push!(braced_tokens, (-1, Int32(-1), semicolon))
    return braced_tokens
end

function extract_parensed(tokens, stack, start, qasm)
    opener  = popfirst!(tokens)
    opener[end] == lparen || throw(QasmParseError("parentethical expression does not open with (", stack, start, qasm))
    openers_met  = 1
    closers_met  = 0
    interior_tokens = Tuple{Int64, Int32, Token}[]
    while closers_met < openers_met && !isempty(tokens)
        next_token      = popfirst!(tokens)
        next_token[end] == lparen && (openers_met += 1)
        next_token[end] == rparen && (closers_met += 1)
        push!(interior_tokens, next_token)
    end
    pop!(interior_tokens) # closing paren
    return interior_tokens
end

function parse_bracketed_expression(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    interior_tokens = extract_braced_block(tokens, stack, start, qasm)
    push!(interior_tokens, (-1, Int32(-1), semicolon))
    return parse_expression(interior_tokens, stack, start, qasm)
end

function parse_paren_expression(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    interior_tokens = extract_parensed(tokens, stack, start, qasm)
    push!(interior_tokens, (-1, Int32(-1), semicolon))
    return parse_expression(interior_tokens, stack, start, qasm)
end

function parse_list_expression(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    expr_list = QasmExpression[]
    while !isempty(tokens) && first(tokens)[end] != semicolon
        tokens[1][end] == comma && popfirst!(tokens)
        next_expr = parse_expression(tokens, stack, start, qasm)
        push!(expr_list, next_expr)
    end
    if length(expr_list) == 1
        return only(expr_list)
    else
        return QasmExpression(:array_literal, expr_list)
    end
end

function parse_literal(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    tokens[1][end] == duration_literal && return parse_duration_literal(popfirst!(tokens), qasm)
    tokens[1][end] == string_token     && return parse_string_literal(popfirst!(tokens), qasm)
    tokens[1][end] == hex              && return parse_hex_literal(popfirst!(tokens), qasm)
    tokens[1][end] == oct              && return parse_oct_literal(popfirst!(tokens), qasm)
    tokens[1][end] == bin              && return parse_bin_literal(popfirst!(tokens), qasm)
    tokens[1][end] == irrational       && return parse_irrational_literal(popfirst!(tokens), qasm)
    tokens[1][end] == boolean          && return parse_boolean_literal(popfirst!(tokens), qasm)
    tokens[1][end] == integer_token    && length(tokens) == 1 && return parse_integer_literal(popfirst!(tokens), qasm)
    tokens[1][end] == float_token      && length(tokens) == 1 && return parse_float_literal(popfirst!(tokens), qasm)
    
    is_float     = tokens[1][end] == float_token
    is_complex   = false
    is_operator  = tokens[2][end] == operator
    is_plusminus = is_operator && parse_identifier(tokens[2], qasm).args[1] ∈ ("+","-")
    is_terminal  = (tokens[2][end] == semicolon || tokens[2][end] == comma || (is_operator && !is_plusminus))
    tokens[1][end] == integer_token && is_terminal && return parse_integer_literal(popfirst!(tokens), qasm)
    tokens[1][end] == float_token   && is_terminal && return parse_float_literal(popfirst!(tokens), qasm)
    splice_end = 1
    if tokens[2][end] == im_token
        is_complex = true
        splice_end = 2
    elseif is_plusminus && tokens[3][end] ∈ (integer_token, float_token) && tokens[4][end] == im_token
        is_complex = true
        is_float |= tokens[3][end] == float_token
        splice_end = 4
    elseif tokens[2][end] ∈ (integer_token, float_token) && tokens[3][end] == im_token # may have absorbed +/- sign
        is_complex = true
        is_float |= tokens[2][end] == float_token
        splice_end = 3
    end
    literal_tokens = splice!(tokens, 1:splice_end)
    raw_literal_string = qasm[literal_tokens[1][1]:literal_tokens[end][1]+literal_tokens[end][2]-1]
    raw_literal = if is_float
            parse(ComplexF64, raw_literal_string)
        elseif is_complex
            parse(Complex{Int}, raw_literal_string)
        else
            parse(Int, raw_literal_string)
        end
    if is_complex # complex float
        return QasmExpression(:complex_literal, raw_literal)
    elseif is_float
        return QasmExpression(:float_literal, raw_literal)
    else
        return QasmExpression(:integer_literal, raw_literal)
    end
    throw(QasmParseError("unable to parse literal", stack, start, qasm)) 
end

function parse_qubit_declaration(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    next_token = tokens[1]
    if next_token[end] == lbracket
        size_tokens = extract_braced_block(tokens, stack, start, qasm)
        size = parse_expression(push!(size_tokens, (-1,-1,semicolon)), stack, start, qasm)
    else
        size = QasmExpression(:integer_literal, 1)
    end
    qubit_name = parse_identifier(popfirst!(tokens), qasm)
    size.args[1] == -1 && (size.args[1] = 1)
    return QasmExpression(:qubit_declaration, qubit_name, size)
end

function parse_gate_mods(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    mod_type = popfirst!(tokens)
    expr = if mod_type[end] == control_mod
               QasmExpression(:control_mod)
           elseif mod_type[end] == negctrl_mod
               QasmExpression(:negctrl_mod)
           elseif mod_type[end] == inverse_mod
               QasmExpression(:inverse_mod)
           elseif mod_type[end] == power_mod
               QasmExpression(:power_mod)
           else
               throw(QasmParseError("cannot parse token of type $(mod_type[end]) as a gate modifier", stack, start, qasm))
           end
    next_token = first(tokens)
    if next_token[end] == lparen
        arg = parse_paren_expression(tokens, stack, start, qasm)
        push!(expr, arg)
        next_token = first(tokens)
    end
    if next_token[end] == at
        popfirst!(tokens)
        next_token = first(tokens)
        if next_token[end] == identifier || next_token[end] == builtin_gate
            push!(expr, parse_expression(tokens, stack, start, qasm))
            return expr
        else
            next_mod_expr = parse_gate_mods(tokens, stack, start, qasm)
            push!(expr, next_mod_expr)
            return expr
        end
    end
end

function parse_expression(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    start_token = popfirst!(tokens)
    next_token  = first(tokens)
    token_name  = QasmExpression(:empty)
    if start_token[end] != classical_type && next_token[end] == lbracket
        name     = parse_identifier(start_token, qasm)
        indices  = parse_expression(tokens, stack, start, qasm)
        expr_indices = length(indices) == 1 ? only(indices) : indices
        token_name   = QasmExpression(:indexed_identifier, name, expr_indices)
    elseif start_token[end] == identifier || start_token[end] == builtin_gate
        token_name = parse_identifier(start_token, qasm)
    elseif start_token[end] == hw_qubit
        token_name = parse_hw_qubit(start_token, qasm)
    elseif start_token[end] == qubit
        token_name = parse_qubit_declaration(tokens, stack, start, qasm)
    elseif start_token[end] == operator
        token_name = parse_identifier(start_token, qasm)
    elseif start_token[end] == break_token
        token_name = QasmExpression(:break)
    elseif start_token[end] == continue_token
        token_name = QasmExpression(:continue)
    elseif start_token[end] == lparen # could be some kind of expression
        token_name = parse_paren_expression(pushfirst!(tokens, start_token), stack, start, qasm) 
    elseif start_token[end] == lbrace # set expression 
        token_name = parse_set_expression(pushfirst!(tokens, start_token), stack, start, qasm) 
    elseif start_token[end] == lbracket
        token_name = parse_bracketed_expression(pushfirst!(tokens, start_token), stack, start, qasm) 
    elseif start_token[end] == classical_type 
        token_name = parse_classical_type(pushfirst!(tokens, start_token), stack, start, qasm)
    elseif start_token[end] ∈ (string_token, integer_token, float_token, hex, oct, bin, irrational, dot, boolean, duration_literal)
        token_name = parse_literal(pushfirst!(tokens, start_token), stack, start, qasm)
    elseif start_token[end] ∈ (mutable, readonly, const_token)
        token_name = parse_identifier(start_token, qasm)
    elseif start_token[end] == dim_token
        raw_dim = qasm[start_token[1]:start_token[1]+start_token[2]-1]
        dim     = replace(replace(raw_dim, " "=>""), "#dim="=>"")
        token_name = QasmExpression(:n_dims, QasmExpression(:integer_literal, parse(Int, dim)))
    end
    head(token_name) == :empty && throw(QasmParseError("unable to parse line with start token $(start_token[end])", stack, start, qasm))
    next_token = first(tokens)
    if next_token[end] == semicolon || next_token[end] == comma || start_token[end] ∈ (lbracket, lbrace)
        expr = token_name
    elseif start_token[end] == integer_token && next_token[end] == irrational # this is banned! 2π is not supported, 2*π is.
        integer_lit = parse_integer_literal(start_token, qasm).args[1]
        irrational_lit = parse_irrational_literal(next_token, qasm).args[1]
        throw(QasmParseError("expressions of form $(integer_lit)$(irrational_lit) are not supported -- you must separate the terms with a * operator.", stack, start, qasm))
    elseif start_token[end] == operator
        unary_op_symbol::Symbol = Symbol(token_name.args[1]::String)
        unary_op_symbol ∈ (:~, :!, :-) || throw(QasmParseError("invalid unary operator $unary_op_symbol.", stack, start, qasm))
        next_expr = parse_expression(tokens, stack, start, qasm)
        # apply unary op to next_expr
        if head(next_expr) ∈ (:identifier, :indexed_identifier, :integer_literal, :float_literal, :string_literal, :irrational_literal, :boolean_literal, :complex_literal, :function_call, :cast, :duration_literal)
            expr = QasmExpression(:unary_op, unary_op_symbol, next_expr)
        elseif head(next_expr) == :binary_op
            # replace first argument
            left_hand_side = next_expr.args[2]::QasmExpression
            new_left_hand_side = QasmExpression(:unary_op, unary_op_symbol, left_hand_side)
            next_expr.args[2] = new_left_hand_side
            expr = next_expr
        end
    elseif next_token[end] == colon
        start = token_name
        popfirst!(tokens)
        second_colon = findfirst(triplet->triplet[end] == colon, tokens)
        if !isnothing(second_colon)
            step_tokens = push!(splice!(tokens, 1:second_colon-1), (-1, Int32(-1), semicolon))
            popfirst!(tokens) # colon
            step = parse_expression(step_tokens, stack, start, qasm)::QasmExpression
        else
            step = QasmExpression(:integer_literal, 1)
        end
        if isempty(tokens) || first(tokens)[end] == semicolon #missing stop
            stop = QasmExpression(:integer_literal, -1)
        else
            stop = parse_expression(tokens, stack, start, qasm)::QasmExpression
        end
        expr = QasmExpression(:range, QasmExpression[start, step, stop])
    elseif next_token[end] == classical_type && start_token[end] ∈ (mutable, readonly, const_token)
        type = parse_classical_type(tokens, stack, start, qasm)
        is_mutable = (start_token[end] == mutable)
        header = is_mutable ? :classical_declaration : :const_declaration
        expr   = QasmExpression(header, type)
        push!(expr, parse_expression(tokens, stack, start, qasm))
    elseif start_token[end] == classical_type && (next_token[end] ∈ (lbracket, identifier))
        expr = QasmExpression(:classical_declaration, token_name)
        push!(expr, parse_expression(tokens, stack, start, qasm))
    elseif start_token[end] == classical_type && next_token[end] == lparen
        expr = QasmExpression(:cast, token_name)
        interior_tokens = extract_parensed(tokens, stack, start, qasm)
        push!(interior_tokens, (-1, Int32(-1), semicolon))
        interior = parse_expression(interior_tokens, stack, start, qasm)
        push!(expr, interior)
    elseif next_token[end] == assignment
        op_token = popfirst!(tokens)
        next_token = first(tokens)
        if next_token[end] ∈ (lparen, lbracket, lbrace, string_token, integer_token, float_token, hex, oct, bin)
            right_hand_side = parse_expression(tokens, stack, start, qasm)
        elseif next_token[end] == measure
            popfirst!(tokens)
            right_hand_side = QasmExpression(:measure, parse_expression(tokens, stack, start, qasm))
        elseif next_token[end] == operator
            unary_op_token = parse_identifier(popfirst!(tokens), qasm)
            next_token = first(tokens)
            unary_right_hand_side = next_token[end] == lparen ? parse_paren_expression(tokens, stack, start, qasm) : parse_expression(tokens, stack, start, qasm)
            right_hand_side = QasmExpression(:unary_op, Symbol(unary_op_token.args[1]::String), unary_right_hand_side)
        else
            right_hand_side = parse_expression(tokens, stack, start, qasm)::QasmExpression
        end
        op_expr = QasmExpression(:binary_op, parse_assignment_op(op_token, qasm), token_name, right_hand_side) 
        expr = QasmExpression(:classical_assignment, op_expr)
    elseif next_token[end] == operator
        op_token = parse_identifier(popfirst!(tokens), qasm)
        right_hand_side = parse_expression(tokens, stack, start, qasm)::QasmExpression
        expr = QasmExpression(:binary_op, Symbol(op_token.args[1]), token_name, right_hand_side)
    else # some kind of function or gate call
        # either a gate call or function call
        arguments  = parse_arguments_list(tokens, stack, start, qasm)
        next_token = first(tokens)
        is_gphase::Bool = (token_name isa QasmExpression && head(token_name) == :identifier && token_name.args[1]::String == "gphase")::Bool
        # this is a gate call with qubit targets
        is_gate_call = next_token[end] == identifier || next_token[end] == hw_qubit || is_gphase
        # this is a function call - unless it is gphase!
        if (next_token[end] == semicolon && !is_gphase)
            popfirst!(tokens)
            expr = QasmExpression(:function_call, token_name, arguments)
        elseif next_token[end] == operator # actually a binary op!
            op_token = parse_identifier(popfirst!(tokens), qasm)
            left_hand_side = QasmExpression(:function_call, token_name, arguments)
            right_hand_side = parse_expression(tokens, stack, start, qasm)
            expr = QasmExpression(:binary_op, Symbol(op_token.args[1]), left_hand_side, right_hand_side)
        else # it's a gate call or gphase
            target_expr = QasmExpression(:qubit_targets, parse_list_expression(tokens, stack, start, qasm))
            expr = QasmExpression(:gate_call, token_name, arguments)
            push!(expr, target_expr)
        end
    end
    return expr
end

function parse_matrix(tokens::Vector{Tuple{Int64, Int32, Token}}, stack, start, qasm)
    inner = extract_braced_block(tokens, stack, start, qasm)
    n_rows = count(triplet->triplet[end] == lbracket, inner)
    matrix = Matrix{QasmExpression}(undef, n_rows, n_rows)
    row = 1
    while !isempty(inner)
        row_tokens = extract_braced_block(inner, stack, start, qasm)
        push!(row_tokens, (-1, Int32(-1), semicolon))
        col = 1
        while !isempty(row_tokens)
            matrix[row, col] = parse_expression(row_tokens, stack, start, qasm)
            col += 1
            next_token = first(row_tokens)
            next_token[end] == comma && popfirst!(row_tokens)
            next_token[end] == semicolon && break
        end
        row += 1
        next_token = first(inner)
        next_token[end] == comma && popfirst!(inner)
        next_token[end] == semicolon && break
    end
    return matrix
end

function parse_qasm(clean_tokens::Vector{Tuple{Int64, Int32, Token}}, qasm::String, root=QasmExpression(:program))
    stack = Stack{QasmExpression}()
    push!(stack, root)
    while !isempty(clean_tokens)
        start, len, token = popfirst!(clean_tokens)
        if token == newline
            continue
        elseif token == version
            closing = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for OPENQASM", stack, start, qasm))
            closing == 1 && throw(QasmParseError("missing version number", stack, start, qasm))
            version_val = parse_literal([popfirst!(clean_tokens)], stack, start, qasm)
            isinteger(version_val.args[1]) || throw(QasmParseError("version number must be an integer", stack, version_start, qasm))
            expr = QasmExpression(:version, QasmExpression(:float_literal, version_val.args[1]))
            push!(stack, expr)
        elseif token == pragma
            closing = findfirst(triplet->triplet[end] == newline, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final newline for #pragma", stack, start, qasm))
            pragma_tokens = splice!(clean_tokens, 1:closing-1)
            push!(stack, parse_pragma(pragma_tokens, stack, start, qasm))
        elseif token == include_token
            closing       = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for include", stack, start, qasm))
            file_name     = popfirst!(clean_tokens)
            popfirst!(clean_tokens) #semicolon
            file_name[end] == string_token || throw(QasmParseError("included filename must be passed as a string", stack, start, qasm))
            file_name_str = replace(qasm[file_name[1]:file_name[1]+file_name[2]-1], "\""=>"", "'"=>"")
            file_contents = read(file_name_str, String)
            file_expr     = parse_qasm(file_contents)
            file_exprs    = file_expr.args
            foreach(ex->push!(stack, ex), file_exprs)
        elseif token == const_token
            closing   = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for const declaration", stack, start, qasm))
            raw_expr = parse_expression(splice!(clean_tokens, 1:closing), stack, start, qasm)
            expr = QasmExpression(:const_declaration, raw_expr.args)
            push!(stack, expr)
        elseif token == classical_type 
            closing   = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for classical declaration", stack, start, qasm))
            line_tokens = pushfirst!(splice!(clean_tokens, 1:closing), (start, len, token))
            expr        = parse_expression(line_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == input
            closing   = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for input", stack, start, qasm))
            input_var = parse_classical_var(splice!(clean_tokens, 1:closing), stack, start, qasm)
            expr      = QasmExpression(:input, input_var...)
            push!(stack, expr)
            popfirst!(clean_tokens) #semicolon
        elseif token == output 
            closing    = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for output", stack, start, qasm))
            output_var = parse_classical_var(splice!(clean_tokens, 1:closing), stack, start, qasm)
            expr       = QasmExpression(:output, output_var...)
            push!(stack, expr)
            popfirst!(clean_tokens) #semicolon
        elseif token == qubit
            closing = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            isnothing(closing) && throw(QasmParseError("missing final semicolon for qubit", stack, start, qasm))
            qubit_tokens = splice!(clean_tokens, 1:closing-1)
            popfirst!(clean_tokens) # semicolon
            expr = parse_qubit_declaration(qubit_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == gate_def
            expr = parse_gate_def(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == function_def
            expr = parse_function_def(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == for_block
            loop_in   = findfirst(triplet->triplet[end] == in_token, clean_tokens)
            isnothing(loop_in) && throw(QasmParseError("for loop variable must have in declaration", stack, start, qasm))
            loop_var  = parse_classical_var(splice!(clean_tokens, 1:loop_in-1), stack, start, qasm)
            popfirst!(clean_tokens) # in
            loop_vals = parse_expression(clean_tokens, stack, start, qasm)
            expr      = parse_for_loop(clean_tokens, loop_var[1], loop_var[2], loop_vals, stack, start, qasm)
            push!(stack, expr)
        elseif token == while_block
            expr = parse_while_loop(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == if_block
            expr = parse_if_block(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == break_token
            push!(stack, QasmExpression(:break))
        elseif token == continue_token
            push!(stack, QasmExpression(:continue))
        elseif token == switch_block
            expr = parse_switch_block(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == line_comment
            eol = findfirst(triplet->triplet[end] == newline, clean_tokens)
            splice!(clean_tokens, 1:eol)
        elseif token == measure 
            eol = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            measure_tokens = splice!(clean_tokens, 1:eol)
            arrow_location = findfirst(triplet->triplet[end] == arrow_token, measure_tokens)
            if !isnothing(arrow_location) # assignment
                targets_tokens = splice!(measure_tokens, 1:arrow_location - 1)
                popfirst!(measure_tokens) # arrow
                targets = parse_expression(push!(targets_tokens, (-1, Int32(-1), semicolon)), stack, start, qasm)
                left_hand_side = parse_expression(measure_tokens, stack, start, qasm)
                right_hand_side = QasmExpression(:measure, targets)
                op_expression = QasmExpression(:binary_op, Symbol("="), left_hand_side, right_hand_side)
                push!(stack, QasmExpression(:classical_assignment, op_expression))
            else
                targets = parse_expression(measure_tokens, stack, start, qasm)
                push!(stack, QasmExpression(:measure, targets))
            end
        elseif token ∈ (negctrl_mod, control_mod, inverse_mod, power_mod)
            gate_mod_tokens = pushfirst!(clean_tokens, (start, len, token))
            expr = parse_gate_mods(gate_mod_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == return_token
            eol = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            return_line_tokens = splice!(clean_tokens, 1:eol)
            line_body  = parse_qasm(return_line_tokens, qasm, QasmExpression(:return))
            line_exprs = collect(Iterators.reverse(line_body))[2:end]
            push!(stack, QasmExpression(:return, line_exprs))
        elseif token == box
            @warn "box expression encountered -- currently `box` is a no-op"
            box_expr = QasmExpression(:box)
            parse_block_body(box_expr, clean_tokens, stack, start, qasm)
            push!(stack, box_expr)
        elseif token == reset_token
            @warn "reset expression encountered -- currently `reset` is a no-op"
            eol = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            reset_tokens = splice!(clean_tokens, 1:eol)
            targets = parse_list_expression(reset_tokens, stack, start, qasm)
            push!(stack, QasmExpression(:reset, targets))
        elseif token == barrier_token
            @warn "barrier expression encountered -- currently `barrier` is a no-op"
            eol = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            barrier_tokens = splice!(clean_tokens, 1:eol)
            targets = parse_list_expression(barrier_tokens, stack, start, qasm)
            push!(stack, QasmExpression(:barrier, targets))
        elseif token == delay_token
            @warn "delay expression encountered -- currently `delay` is a no-op"
            eol = findfirst(triplet->triplet[end] == semicolon, clean_tokens)
            delay_tokens = splice!(clean_tokens, 1:eol)
            delay_expr   = QasmExpression(:delay)
            # format is delay[duration]; or delay[duration] targets;
            delay_duration = extract_braced_block(delay_tokens, stack, start, qasm)
            push!(delay_expr, QasmExpression(:duration, parse_expression(delay_duration, stack, start, qasm)))
            target_expr = QasmExpression(:targets)
            if first(delay_tokens)[end] != semicolon # targets present
                targets = parse_list_expression(delay_tokens, stack, start, qasm)
                push!(target_expr, targets)
            end
            push!(delay_expr, target_expr)
            push!(stack, delay_expr)
        elseif token == end_token
            push!(stack, QasmExpression(:end))
        elseif token == identifier || token == builtin_gate
            clean_tokens = pushfirst!(clean_tokens, (start, len, token))
            expr = parse_expression(clean_tokens, stack, start, qasm)
            push!(stack, expr)
        elseif token == forbidden_keyword
            token_id = name(parse_identifier((start, len, token), qasm))
            throw(QasmParseError("keyword $token_id not supported.", stack, start, qasm))
        end
    end
    return stack
end
function parse_qasm(qasm::String, root=QasmExpression(:program))::QasmExpression
    raw_tokens   = tokenize(Token, qasm)
    clean_tokens = filter(triplet->triplet[3] ∉ (spaces, block_comment), collect(raw_tokens))
    # add a final newline in case one is missing 
    clean_tokens[end][end] == newline || push!(clean_tokens, (-1, Int32(-1), newline))
    stack = parse_qasm(clean_tokens, qasm, root)
    stack_exprs = convert(Vector{QasmExpression}, collect(Iterators.reverse(stack)))::Vector{QasmExpression}
    foreach(ex->push!(stack_exprs[1], ex), stack_exprs[2:end])
    return stack_exprs[1] 
end

mutable struct ClassicalVariable
    name::String
    type
    val
    is_const::Bool
end

struct Qubit 
    name::String
    size::Int
end

const InstructionArgument = Union{Symbol, Dates.Period, Real, Matrix{ComplexF64}}
const CircuitInstruction = @NamedTuple begin type::String; arguments::Vector{InstructionArgument}; targets::Vector{Int}; controls::Vector{Pair{Int, Int}}; exponent::Float64 end
const CircuitResult = @NamedTuple begin type::Symbol; operator::Vector{Union{String, Matrix{ComplexF64}}}; targets::Vector{Int}; states::Vector{String}; end

abstract type AbstractGateDefinition end

struct GateDefinition <: AbstractGateDefinition
    name::String
    arguments::Vector{String}
    qubit_targets::Vector{String} # keep this as string to support splatting
    body::QasmExpression
end

struct BuiltinGateDefinition <: AbstractGateDefinition
    name::String
    arguments::Vector{String}
    qubit_targets::Vector{String} # keep this as string to support splatting
    body::CircuitInstruction
end

struct FunctionDefinition
    name::String
    arguments::QasmExpression
    body::Vector{QasmExpression}
    return_type
end
FunctionDefinition(name::String, arguments::QasmExpression, body::QasmExpression, return_type) = FunctionDefinition(name, arguments, [body], return_type)

struct QasmVisitorError <: Exception
    message::String
    alternate_type::String
end
QasmVisitorError(message::String) = QasmVisitorError(message, "")
function Base.showerror(io::IO, err::QasmVisitorError)
    print(io, "QasmVisitorError: ")
    print(io, err.message)
end

abstract type AbstractVisitor end

include("builtin_gates.jl")

mutable struct QasmProgramVisitor <: AbstractVisitor
    inputs::Dict{String, Any}
    classical_defs::Dict{String, ClassicalVariable}
    function_defs::Dict{String, FunctionDefinition}
    gate_defs::Dict{String, AbstractGateDefinition}
    qubit_defs::Dict{String, Qubit}
    qubit_mapping::Dict{String, Vector{Int}}
    qubit_count::Int
    instructions::Vector{CircuitInstruction}
    results::Vector{CircuitResult}
    function QasmProgramVisitor(inputs::Dict{String, <:Any} = Dict{String, Any}())
        new(inputs,
            Dict{String, ClassicalVariable}(),
            Dict{String, FunctionDefinition}(),
            builtin_gates(),
            Dict{String, Qubit}(),
            Dict{String, Vector{Int}}(),
            0,
            CircuitInstruction[],
            CircuitResult[],
           )
    end
end

mutable struct QasmGateDefVisitor <: AbstractVisitor
    parent::AbstractVisitor
    classical_defs::Dict{String, ClassicalVariable}
    qubit_defs::Dict{String, Qubit}
    qubit_mapping::Dict{String, Vector{Int}}
    qubit_count::Int
    instructions::Vector{CircuitInstruction}
    function QasmGateDefVisitor(parent::AbstractVisitor, declared_arguments::Vector{String}, provided_arguments::QasmExpression, gate_qubits::Vector{String})
        qubit_defs    = Dict(q=>Qubit(q, 1) for q in gate_qubits)
        qubit_mapping = Dict(gate_qubits[ix+1]=>[ix] for ix in 0:length(gate_qubits)-1)
        for ix in 0:length(gate_qubits)-1
            qubit_mapping[gate_qubits[ix+1] * "[0]"] = [ix]
        end
        v = new(parent,
                deepcopy(classical_defs(parent)),
                qubit_defs,
                qubit_mapping,
                length(gate_qubits),
                CircuitInstruction[],
               )
        for (arg_name, arg_value) in zip(declared_arguments, parent(provided_arguments))
            classical_defs(v)[arg_name] = ClassicalVariable(arg_name, Real, arg_value, true)
        end
        return v
    end
end

mutable struct QasmForLoopVisitor <: AbstractVisitor
    parent::AbstractVisitor
    classical_defs::Dict{String, ClassicalVariable}
    QasmForLoopVisitor(parent::AbstractVisitor) = new(parent, classical_defs(parent))
end

mutable struct QasmWhileLoopVisitor <: AbstractVisitor
    parent::AbstractVisitor
    QasmWhileLoopVisitor(parent::AbstractVisitor) = new(parent)
end

mutable struct QasmFunctionVisitor <: AbstractVisitor
    parent::AbstractVisitor
    classical_defs::Dict{String, ClassicalVariable}
    qubit_defs::Dict{String, Qubit}
    qubit_mapping::Dict{String, Vector{Int}}
    qubit_count::Int
    instructions::Vector{CircuitInstruction}
    function QasmFunctionVisitor(parent::AbstractVisitor, declared_arguments::Vector{QasmExpression}, provided_arguments::Vector{QasmExpression})
        v = new(parent, 
            classical_defs(parent),
            deepcopy(parent.qubit_defs),
            deepcopy(parent.qubit_mapping),
            qubit_count(parent),
            CircuitInstruction[],
           )
        arg_map = Dict(zip(declared_arguments, provided_arguments))
        for arg in declared_arguments
            if head(arg) ∈ (:const_declaration, :classical_declaration)
                new_val = parent(arg_map[arg])
                if head(arg.args[2]) != :classical_assignment
                    arg_id = pop!(arg)
                    push!(arg, QasmExpression(:classical_assignment, QasmExpression(:binary_op, Symbol("="), arg_id, new_val)))
                else
                    arg.args[2].args[1].args[end] = new_val
                end
            end
            v(arg)
        end
        return v
    end
end
function QasmFunctionVisitor(parent::AbstractVisitor, declared_arguments::Vector{QasmExpression}, provided_arguments::QasmExpression)
    head(provided_arguments) == :array_literal && return QasmFunctionVisitor(parent, declared_arguments, convert(Vector{QasmExpression}, provided_arguments.args))
    QasmFunctionVisitor(parent, declared_arguments, [provided_arguments])
end
function QasmFunctionVisitor(parent::AbstractVisitor, declared_arguments::QasmExpression, provided_arguments)
    head(declared_arguments) == :array_literal && return QasmFunctionVisitor(parent, convert(Vector{QasmExpression}, declared_arguments.args), provided_arguments)
    QasmFunctionVisitor(parent, [declared_arguments], provided_arguments)
end
Base.parent(v::AbstractVisitor) = v.parent

hasgate(v::AbstractVisitor, gate_name::String)    = hasgate(parent(v), gate_name)
hasgate(v::QasmProgramVisitor, gate_name::String) = haskey(v.gate_defs, gate_name)
gate_defs(v::AbstractVisitor)    = gate_defs(parent(v))
gate_defs(v::QasmProgramVisitor) = v.gate_defs

function_defs(v::QasmProgramVisitor) = v.function_defs
function_defs(v::AbstractVisitor)    = function_defs(parent(v))

hasfunction(v::AbstractVisitor, function_name::String)    = haskey(function_defs(v), function_name)

qubit_defs(v::AbstractVisitor)     = qubit_defs(parent(v))
qubit_defs(v::QasmFunctionVisitor) = v.qubit_defs
qubit_defs(v::QasmProgramVisitor)  = v.qubit_defs

qubit_mapping(v::AbstractVisitor)     = qubit_mapping(parent(v))
qubit_mapping(v::QasmProgramVisitor)  = v.qubit_mapping
qubit_mapping(v::QasmFunctionVisitor) = v.qubit_mapping
qubit_mapping(v::QasmGateDefVisitor)  = v.qubit_mapping

qubit_count(o::String)              = length(o)
qubit_count(o::Matrix)              = Int(log2(size(o, 1)))
qubit_count(v::AbstractVisitor)     = qubit_count(parent(v))
qubit_count(v::QasmProgramVisitor)  = v.qubit_count
qubit_count(v::QasmFunctionVisitor) = v.qubit_count
qubit_count(v::QasmGateDefVisitor)  = v.qubit_count

classical_defs(v::AbstractVisitor)     = classical_defs(parent(v))
classical_defs(v::QasmProgramVisitor)  = v.classical_defs
classical_defs(v::QasmGateDefVisitor)  = v.classical_defs
classical_defs(v::QasmFunctionVisitor) = v.classical_defs

instructions(v::AbstractVisitor)     = instructions(parent(v))
instructions(v::QasmProgramVisitor)  = v.instructions
instructions(v::QasmGateDefVisitor)  = v.instructions
instructions(v::QasmFunctionVisitor) = v.instructions

results(v::AbstractVisitor) = results(parent(v))
results(v::QasmProgramVisitor) = v.results

Base.push!(v::AbstractVisitor, ixs::Vector{CircuitInstruction}) = append!(instructions(v), ixs)
Base.push!(v::AbstractVisitor, ix::CircuitInstruction)          = push!(instructions(v), ix)

Base.push!(v::AbstractVisitor, rts::Vector{CircuitResult}) = append!(results(v), rts)
Base.push!(v::AbstractVisitor, rt::CircuitResult)          = push!(results(v), rt)

function evaluate_unary_op(op::Symbol, arg)
    op == :! && return !arg
    op == :~ && return .!arg
    op == :- && return -arg 
end
function evaluate_unary_op(op::Symbol, arg::BitVector)
    op == :! && return !any(arg)
    op == :~ && return .!arg
end

# semgrep rules can't handle this macro properly yet
# nosemgrep
function evaluate_binary_op(op::Symbol, @nospecialize(lhs), @nospecialize(rhs))
    op == :< && return lhs < rhs
    op == :> && return lhs > rhs
    op == :<= && return lhs <= rhs
    op == :>= && return lhs >= rhs
    op == Symbol("=") && return rhs
    op == Symbol("!=") && return lhs != rhs
    op == Symbol("==") && return lhs == rhs
    op == :+ && return lhs + rhs
    op == :- && return lhs - rhs
    op == :* && return lhs * rhs
    op == :/ && return lhs / rhs
    op == :% && return lhs % rhs
    op == Symbol("<<") && return lhs << rhs
    op == Symbol(">>") && return lhs >> rhs
    op == Symbol("**") && return lhs ^ rhs
    op == Symbol("&&") && return lhs && rhs
    op == Symbol("||") && return lhs || rhs
    op == :| && return lhs .| rhs
    op == :& && return lhs .& rhs
    op == :^ && return lhs .⊻ rhs
end

function name(expr::QasmExpression)::String
    head(expr) == :identifier         && return expr.args[1]::String
    head(expr) == :indexed_identifier && return name(expr.args[1]::QasmExpression)
    head(expr) == :qubit_declaration  && return name(expr.args[1]::QasmExpression)
    head(expr) == :classical_declaration && return name(expr.args[2]::QasmExpression)
    head(expr) == :input              && return name(expr.args[2]::QasmExpression)
    head(expr) == :function_call      && return name(expr.args[1]::QasmExpression)
    head(expr) == :gate_call          && return name(expr.args[1]::QasmExpression)
    head(expr) == :gate_definition    && return name(expr.args[1]::QasmExpression)
    head(expr) == :classical_assignment && return name(expr.args[1].args[2]::QasmExpression)
    head(expr) == :hw_qubit           && return replace(expr.args[1], "\$"=>"")
    throw(QasmVisitorError("name not defined for expressions of type $(head(expr))"))
end

function evaluate_modifiers(v::V, expr::QasmExpression) where {V<:AbstractVisitor}
    if head(expr) == :power_mod
        pow_expr = QasmExpression(:pow, v(expr.args[1]::QasmExpression))
        return (pow_expr, expr.args[2])
    elseif head(expr) == :inverse_mod
        return (QasmExpression(:inv), expr.args[1]::QasmExpression)
    elseif head(expr) ∈ (:control_mod, :negctrl_mod)
        has_argument = length(expr.args) > 1
        if has_argument
            arg_val::Int = v(first(expr.args)::QasmExpression)::Int
            isinteger(arg_val) || throw(QasmVisitorError("cannot apply non-integer ($arg_val) number of controls or negcontrols."))
            true_inner = expr.args[2]::QasmExpression
            inner = QasmExpression(head(expr), true_inner)
            while arg_val > 2
                inner = QasmExpression(head(expr), inner)
                arg_val -= 1
            end
        else
            inner = expr.args[1]::QasmExpression
        end
        new_head = head(expr) == :control_mod ? :ctrl : :negctrl
        return (QasmExpression(new_head), inner)
    end
end

# nosemgrep
function _evaluate_qubits(::Val{:identifier}, v, qubit_expr::QasmExpression)::Vector{Int}
    qubit_name = name(qubit_expr)
    mapping    = qubit_mapping(v)::Dict{String, Vector{Int}}
    haskey(mapping, qubit_name) || throw(QasmVisitorError("Missing input variable '$qubit_name'.", "NameError"))
    return mapping[qubit_name]
end

# nosemgrep
function _evaluate_qubits(::Val{:indexed_identifier}, v, qubit_expr::QasmExpression)::Vector{Int}
    qubit_name = name(qubit_expr)
    mapping    = qubit_mapping(v)::Dict{String, Vector{Int}}
    haskey(mapping, qubit_name) || throw(QasmVisitorError("Missing input variable '$qubit_name'.", "NameError"))
    qubit_ix   = v(qubit_expr.args[2]::QasmExpression)
    qubits     = Iterators.flatmap(qubit_ix) do rq
        haskey(mapping, qubit_name * "[$rq]") || throw(QasmVisitorError("Invalid qubit index '$rq' in '$qubit_name'.", "IndexError"))
        return mapping[qubit_name * "[$rq]"]
    end
    return collect(qubits)
end
_evaluate_qubits(::Val{:array_literal}, v, qubit_expr::QasmExpression)::Vector{Int} = collect(Iterators.flatmap(expr->_evaluate_qubits(Val(head(expr)), v, expr), qubit_expr.args))
_evaluate_qubits(::Val{:hw_qubit}, v, qubit_expr::QasmExpression)::Vector{Int} = Int[v(qubit_expr)::Int]
_evaluate_qubits(val, v, qubit_expr) = throw(QasmVisitorError("unable to evaluate qubits for expression $qubit_expr."))

function evaluate_qubits(v::AbstractVisitor, qubit_targets::Vector)::Vector{Int}
    final_qubits = map(qubit_expr->_evaluate_qubits(Val(head(qubit_expr)), v, qubit_expr), qubit_targets)
    return vcat(final_qubits...)
end
evaluate_qubits(v::AbstractVisitor, qubit_targets::QasmExpression) = evaluate_qubits(v::AbstractVisitor, [qubit_targets])

function remap(ix, target_mapper::Dict{Int, Int})
    mapped_targets  = map(t->getindex(target_mapper, t), ix.targets)
    mapped_controls = map(c->getindex(target_mapper, c[1])=>c[2], ix.controls)
    return (type=ix.type, arguments=ix.arguments, targets=mapped_targets, controls=mapped_controls, exponent=ix.exponent)
end
function bind_arguments!(ix::CircuitInstruction, argument_values::Dict{Symbol, <:Real})
    new_arguments = InstructionArgument[get(argument_values, arg, arg) for arg in ix.arguments]
    return (type=ix.type, arguments=new_arguments, targets=ix.targets, controls=ix.controls, exponent=ix.exponent)
end

function process_gate_arguments(v::AbstractVisitor, gate_name::String, defined_arguments::Vector{String}, called_arguments::QasmExpression, gate_body::Vector{CircuitInstruction})
    def_has_arguments  = !isempty(defined_arguments)
    call_has_arguments = !isempty(v(called_arguments))
    if def_has_arguments ⊻ call_has_arguments
        def_has_arguments && throw(QasmVisitorError("gate $gate_name requires arguments but none were provided.")) 
        call_has_arguments && throw(QasmVisitorError("gate $gate_name does not accept arguments but arguments were provided."))
    end
    if def_has_arguments
        evaled_args     = v(called_arguments)
        argument_values = Dict{Symbol, Real}(Symbol(arg_name)=>argument for (arg_name, argument) in zip(defined_arguments, evaled_args))
        return map(ix->bind_arguments!(ix, argument_values), gate_body)
    else
        return deepcopy(gate_body)
    end
end

function handle_gate_modifiers(ixs, mods::Vector{QasmExpression}, control_qubits::Vector{Int}, is_gphase::Bool)
    for mod in Iterators.reverse(mods)
        control_qubit = head(mod) ∈ (:negctrl, :ctrl) ? pop!(control_qubits) : -1
        for (ii, ix) in enumerate(ixs)
            if head(mod) == :pow
                ixs[ii] = (type=ix.type, arguments=ix.arguments, targets=ix.targets, controls=ix.controls, exponent=ix.exponent*mod.args[1])
            elseif head(mod) == :inv
                ixs[ii] = (type=ix.type, arguments=ix.arguments, targets=ix.targets, controls=ix.controls, exponent=-ix.exponent)
            # need to handle "extra" target
            elseif head(mod) ∈ (:negctrl, :ctrl)
                bit = head(mod) == :ctrl ? 1 : 0
                if is_gphase
                    ixs[ii] = (type=ix.type, arguments=ix.arguments, targets=ix.targets, controls=pushfirst!(ix.controls, control_qubit=>bit), exponent=ix.exponent)
                else
                    ixs[ii] = (type=ix.type, arguments=ix.arguments, targets=pushfirst!(ix.targets, control_qubit), controls=pushfirst!(ix.controls, control_qubit=>bit), exponent=ix.exponent)
                end
            end
        end
        head(mod) == :inv && reverse!(ixs)
    end
    return ixs
end

function splat_gate_targets(gate_targets::Vector{Vector{Int}})
    target_lengths::Vector{Int} = Int[length(t) for t in gate_targets]
    longest = maximum(target_lengths)
    must_splat::Bool = any(len->len!=1 || len != longest, target_lengths)
    !must_splat && return longest, gate_targets
    for target_ix in 1:length(gate_targets)
        if target_lengths[target_ix] == 1
            append!(gate_targets[target_ix], fill(only(gate_targets[target_ix]), longest-1))
        end
    end
    return longest, gate_targets
end

function visit_gphase_call(v::AbstractVisitor, program_expr::QasmExpression)
    gate_name::String = program_expr.args[1].args[1]
    is_gphase = gate_name == "gphase"
    call_targets::Vector{QasmExpression}  = program_expr.args[3].args
    provided_args::Vector{QasmExpression} = program_expr.args[2].args
    has_modifiers = length(program_expr.args) == 4
    n_called_with::Int  = qubit_count(v)
    n_defined_with::Int = n_called_with
    gate_targets::Vector{Int} = collect(0:n_called_with-1)
    provided_arg::QasmExpression = only(program_expr.args[2].args)
    evaled_arg        = Float64(v(provided_arg))
    applied_arguments = CircuitInstruction[(type="gphase", arguments=[evaled_arg], targets=gate_targets, controls=Pair{Int,Int}[], exponent=1.0)]
    mods::Vector{QasmExpression} = length(program_expr.args) == 4 ? program_expr.args[4].args : QasmExpression[]
    applied_arguments = handle_gate_modifiers(applied_arguments, mods, deepcopy(gate_targets), true)
    target_mapper     = Dict{Int, Int}(g_ix=>gate_targets[g_ix+1][1] for g_ix in 0:n_called_with-1)
    push!(v, map(ix->remap(ix, target_mapper), applied_arguments))
    return
end

function visit_gate_call(v::AbstractVisitor, program_expr::QasmExpression)
    gate_name        = name(program_expr)::String
    raw_call_targets = program_expr.args[3]::QasmExpression
    call_targets::Vector{QasmExpression}  = convert(Vector{QasmExpression}, head(raw_call_targets.args[1]) == :array_literal ? raw_call_targets.args[1].args : raw_call_targets.args)::Vector{QasmExpression}
    provided_args    = isempty(program_expr.args[2].args) ? QasmExpression(:empty) : only(program_expr.args[2].args)::QasmExpression
    has_modifiers    = length(program_expr.args) == 4
    hasgate(v, gate_name) || throw(QasmVisitorError("gate $gate_name not defined!"))
    gate_def         = gate_defs(v)[gate_name]
    gate_def_v       = QasmGateDefVisitor(v, gate_def.arguments, provided_args, gate_def.qubit_targets)
    gate_def_v(deepcopy(gate_def.body))
    gate_ixs         = instructions(gate_def_v)
    gate_targets     = Vector{Int}[evaluate_qubits(v, call_target)::Vector{Int} for call_target in call_targets]
    n_called_with    = length(gate_targets)
    n_defined_with   = length(gate_def.qubit_targets)
    # cases like `ccnot qs`;
    if n_called_with < n_defined_with && length(gate_targets[1]) == n_defined_with
        n_called_with = length(gate_targets[1])
        gate_targets  = Vector{Int}[[gt] for gt in gate_targets[1]]
    end
    applied_arguments = process_gate_arguments(v, gate_name, gate_def.arguments, provided_args, gate_ixs)
    control_qubits::Vector{Int} = collect(0:(n_called_with-n_defined_with)-1)
    mods::Vector{QasmExpression} = length(program_expr.args) == 4 ? convert(Vector{QasmExpression}, program_expr.args[4].args) : QasmExpression[]
    if !isempty(control_qubits)
        modifier_remap = Dict{Int, Int}(old_qubit=>(old_qubit + length(control_qubits)) for old_qubit in 0:length(gate_def.qubit_targets))
        for ii in 1:length(applied_arguments)
            applied_arguments[ii] = remap(applied_arguments[ii], modifier_remap)
        end
    end
    applied_arguments     = handle_gate_modifiers(applied_arguments, mods, control_qubits, false)
    longest, gate_targets = splat_gate_targets(gate_targets)
    for splatted_ix in 1:longest
        target_mapper = Dict{Int, Int}(g_ix=>gate_targets[g_ix+1][splatted_ix] for g_ix in 0:n_called_with-1)
        push!(v, map(ix->remap(ix, target_mapper), applied_arguments))
    end
    return
end

(v::AbstractVisitor)(i::Number) = i
(v::AbstractVisitor)(i::String) = i
(v::AbstractVisitor)(i::BitVector) = i
(v::AbstractVisitor)(i::NTuple{N,<:Number}) where {N} = i
(v::AbstractVisitor)(i::Vector{<:Number}) = i
(v::AbstractVisitor)(program_exprs::Vector) = map(v, program_exprs)
(v::QasmGateDefVisitor)(ix::CircuitInstruction) = push!(v, ix)
function (v::AbstractVisitor)(program_expr::QasmExpression)
    var_name::String = ""
    if head(program_expr) == :program
        for expr in program_expr.args
            head(expr) == :end && return
            v(expr)
        end
    elseif head(program_expr) == :scope
        for expr in program_expr.args
            head(expr) == :end && return
            head(expr) == :continue && return :continue
            head(expr) == :break && return :break
            v(expr)
        end
    elseif head(program_expr) == :version
        return v
    elseif head(program_expr) == :reset
        targets = program_expr.args[1]::QasmExpression
        target_qubits = evaluate_qubits(v, targets)
        ixs = [(type="reset", arguments=InstructionArgument[], targets=[t], controls=Pair{Int, Int}[], exponent=1.0) for t in target_qubits]
        push!(v, ixs)
        return v
    elseif head(program_expr) == :barrier
        targets = program_expr.args[1]::QasmExpression
        target_qubits = evaluate_qubits(v, targets)
        ixs = [(type="barrier", arguments=InstructionArgument[], targets=[t], controls=Pair{Int, Int}[], exponent=1.0) for t in target_qubits]
        push!(v, ixs)
        return v
    elseif head(program_expr) == :delay
        duration_expr = program_expr.args[1].args[1]::QasmExpression
        targets       = program_expr.args[2].args[1]::QasmExpression
        target_qubits = evaluate_qubits(v, targets)
        duration      = v(duration_expr)
        ixs = [(type="delay", arguments=InstructionArgument[duration], targets=[t], controls=Pair{Int, Int}[], exponent=1.0) for t in target_qubits]
        push!(v, ixs)
        return v
    elseif head(program_expr) == :stretch
        return v
    elseif head(program_expr) == :duration
        return v
    elseif head(program_expr) == :input
        var_name = name(program_expr)
        var_type = program_expr.args[1].args[1]
        haskey(v.inputs, var_name) || throw(QasmVisitorError("Missing input variable '$var_name'.", "NameError"))
        var      = ClassicalVariable(var_name, var_type, v.inputs[var_name], true)
        v.classical_defs[var_name] = var
        return v
    elseif head(program_expr) ∈ (:continue, :break)
        v isa Union{QasmForLoopVisitor, QasmWhileLoopVisitor} && return head(program_expr)
        throw(QasmVisitorError(string(head(program_expr)) * " statement encountered outside a loop."))
    elseif head(program_expr) == :for
        for_v = QasmForLoopVisitor(v)
        for_loop             = convert(Vector{QasmExpression}, program_expr.args)
        loop_variable_type   = for_loop[1].args[1]
        loop_variable_name   = for_loop[2].args[1]::String
        loop_variable_values = for_v(for_loop[3])
        loop_body            = for_loop[4]::QasmExpression
        for loop_value in loop_variable_values
            loop_variable = ClassicalVariable(loop_variable_name, loop_variable_type, loop_value, false)
            for_v.classical_defs[loop_variable_name] = loop_variable
            if head(loop_body) == :scope
                for expr in convert(Vector{QasmExpression}, loop_body.args)
                    rt = for_v(expr)
                    rt == :continue && break
                    if rt == :break
                        delete!(classical_defs(v), loop_variable_name)
                        return v
                    end
                end
            else
                for_v(loop_body)
            end
        end
        delete!(classical_defs(v), loop_variable_name)
    elseif head(program_expr) == :switch
        case_val = v(program_expr.args[1])
        all_cases = convert(Vector{QasmExpression}, program_expr.args[2:end])
        default = findfirst(expr->head(expr) == :default, all_cases)
        case_found = false
        for case in all_cases
            if head(case) == :case && case_val ∈ v(case.args[1])
                case_found = true
                foreach(v, convert(Vector{QasmExpression}, case.args[2:end]))
                break
            end
        end
        if !case_found
            isnothing(default) && throw(QasmVisitorError("no case matched and no default defined."))
            foreach(v, convert(Vector{QasmExpression}, all_cases[default].args))
        end
    elseif head(program_expr) == :identifier
        id_name = name(program_expr)
        haskey(classical_defs(v), id_name) && return classical_defs(v)[id_name].val
        haskey(qubit_mapping(v), id_name) && return evaluate_qubits(v, program_expr)
        throw(QasmVisitorError("no identifier $id_name defined."))
    elseif head(program_expr) == :indexed_identifier
        identifier_name = name(program_expr)
        if haskey(classical_defs(v), identifier_name)
            var = classical_defs(v)[identifier_name]
            ix  = v(program_expr.args[2]::QasmExpression)
            if ix isa StepRange && ix.step > 0 && ix.stop < ix.start # -1 in place of end
                new_stop = var.type isa SizedNumber || var.type isa SizedBitVector ? v(var.type.size) : length(var.val)
                ix = StepRange(ix.start, ix.step, new_stop-1)
            end
            flat_ix = mapreduce(ix_ -> ix_ .+ 1, vcat, ix)
            if var.type isa SizedInt || var.type isa SizedUInt
                n_bits::Int = v(var.type.size)::Int
                int_val     = convert(Int, var.val)::Int
                values      = Int[(int_val >> (n_bits - index)) & 1 for index in flat_ix]
                return length(flat_ix) == 1 ? values[1] : values
            else
                return length(flat_ix) == 1 ? var.val[only(flat_ix)] : var.val[flat_ix]
            end
        elseif haskey(qubit_mapping(v), identifier_name)
            return evaluate_qubits(v, program_expr)
        else
            throw(QasmVisitorError("no identifier $identifier_name defined."))
        end
    elseif head(program_expr) == :if
        condition_value = v(program_expr.args[1]) > 0
        has_else  = findfirst(expr->head(expr) == :else, convert(Vector{QasmExpression}, program_expr.args))
        last_expr = !isnothing(has_else) ? length(program_expr.args) - 1 : length(program_expr.args)
        if condition_value
            for expr in program_expr.args[2:last_expr]
                rt = v(expr)
                rt == :continue && return :continue
                rt == :break && return :break
            end
        elseif !isnothing(has_else)
            for expr in program_expr.args[has_else].args
                rt = v(expr)
                rt == :continue && return :continue
                rt == :break && return :break
            end
        end
    elseif head(program_expr) == :while
        while_v = QasmWhileLoopVisitor(v)
        condition_value = v(program_expr.args[1]) > 0
        loop_body = program_expr.args[2]
        while condition_value
            if head(loop_body) == :scope
                for expr in loop_body.args
                    rt = while_v(expr)
                    rt == :continue && break
                    rt == :break && return v
                end
            else
                while_v(loop_body)
            end
            condition_value = while_v(program_expr.args[1])
        end
    elseif head(program_expr) == :classical_assignment
        op = program_expr.args[1].args[1]::Symbol
        left_hand_side  = program_expr.args[1].args[2]::QasmExpression
        right_hand_side = program_expr.args[1].args[3]
        var_name  = name(left_hand_side)::String
        right_val = v(right_hand_side)
        left_val  = v(left_hand_side)
        classical_defs(v)[var_name].is_const && throw(QasmVisitorError("cannot reassign value of const variable!"))
        if head(left_hand_side) == :identifier
            var = classical_defs(v)[var_name]
            var_type = var.type
            if var_type isa SizedBitVector && right_val isa AbstractString # bitstring literal
                cleaned_val::String = replace(right_val, "\""=>"")
                bit_right = BitVector(tryparse(Int, "$b") for b in cleaned_val)
                new_val = evaluate_binary_op(op, left_val, bit_right)
            else
                new_val = evaluate_binary_op(op, left_val, right_val)
            end
            var.val = new_val 
        elseif head(left_hand_side) == :indexed_identifier
            inds = v(left_hand_side.args[2])
            var  = classical_defs(v)[var_name]
            var_type = var.type
            if inds isa StepRange && inds.step > 0 && inds.stop < inds.start # -1 in place of end
                new_stop = var.type isa SizedNumber ? v(var.type.size) - 1 : length(var.val) - 1 
                inds = StepRange(inds.start, inds.step, new_stop)
            end
            inds = inds .+ 1
            if var_type isa SizedBitVector && right_val isa AbstractString # bitstring literal
                cleaned_val = replace(right_val, "\""=>"")
                bit_right = BitVector(tryparse(Int, "$b") for b in cleaned_val)
                new_val = evaluate_binary_op(op, left_val, bit_right)
            else
                new_val = evaluate_binary_op(op, left_val, right_val)
            end
            if length(inds) > 1
                var.val[inds] .= new_val 
            else
                var.val[inds] = new_val 
            end
        end
    elseif head(program_expr) == :classical_declaration
        var_type = program_expr.args[1].args[1]
        init = if var_type isa SizedNumber
                undef
            elseif var_type isa SizedArray
                fill(undef, v(var_type.size))
            elseif var_type isa SizedBitVector
                falses(max(0, v(var_type.size)))
            end
        # no initial value
        if head(program_expr.args[2]) == :identifier
            var_name = name(program_expr.args[2])
            v.classical_defs[var_name] = ClassicalVariable(var_name, var_type, init, false)
        elseif head(program_expr.args[2]) == :classical_assignment
            op, left_hand_side, right_hand_side = program_expr.args[2].args[1].args
            var_name = name(left_hand_side)
            v.classical_defs[var_name] = ClassicalVariable(var_name, var_type, init, false)
            v(program_expr.args[2])
        end
    elseif head(program_expr) == :const_declaration
        head(program_expr.args[2]) == :classical_assignment || throw(QasmVisitorError("const declaration must assign an initial value."))
        var_type = program_expr.args[1].args[1]
        init = if var_type isa SizedNumber
                undef
            elseif var_type isa SizedArray
                fill(undef, v(var_type.size))
            elseif var_type isa SizedBitVector
                falses(max(0, v(var_type.size)))
            end
        op, left_hand_side, right_hand_side = program_expr.args[2].args[1].args
        var_name = name(left_hand_side)
        v.classical_defs[var_name] = ClassicalVariable(var_name, var_type, init, false)
        v(program_expr.args[2])
        v.classical_defs[var_name] = ClassicalVariable(var_name, var_type, v.classical_defs[var_name].val, true)
    elseif head(program_expr) == :qubit_declaration
        qubit_name::String = name(program_expr)
        qubit_size::Int = v(program_expr.args[2])::Int
        qubit_defs(v)[qubit_name] = Qubit(qubit_name, qubit_size)
        qubit_mapping(v)[qubit_name] = collect(qubit_count(v) : qubit_count(v) + qubit_size - 1)
        for qubit_i in 0:qubit_size-1
            qubit_mapping(v)["$qubit_name[$qubit_i]"] = [qubit_count(v) + qubit_i]
        end
        v.qubit_count += qubit_size
    elseif head(program_expr) ∈ (:power_mod, :inverse_mod, :control_mod, :negctrl_mod)
        mods = QasmExpression(:modifiers)
        mod_expr, inner = evaluate_modifiers(v, program_expr)
        push!(mods, mod_expr)
        while head(inner) != :gate_call # done
            mod_expr, inner = evaluate_modifiers(v, inner)
            push!(mods, mod_expr)
        end
        push!(inner, mods)
        v(inner)
    elseif head(program_expr) == :gate_call
        gate_name = name(program_expr)
        is_gphase = gate_name == "gphase"
        if is_gphase
            visit_gphase_call(v, program_expr)
        else
            visit_gate_call(v, program_expr)
        end
    elseif head(program_expr) == :box
        foreach(v, program_expr.args)
    elseif head(program_expr) == :gate_definition
        gate_def         = program_expr.args
        gate_name        = name(program_expr)
        gate_arguments   = gate_def[2]::QasmExpression
        gate_def_targets = gate_def[3]::QasmExpression
        gate_body        = gate_def[4]::QasmExpression
        single_argument  = !isempty(gate_arguments.args) && head(gate_arguments.args[1]) == :array_literal
        argument_exprs   = single_argument ? gate_arguments.args[1].args::Vector{Any} : gate_arguments.args::Vector{Any}
        argument_names   = String[arg.args[1] for arg::QasmExpression in argument_exprs]
        single_target    = head(gate_def_targets.args[1]) == :array_literal
        qubit_targets    = single_target ? map(name, gate_def_targets.args[1].args)::Vector{String} : map(name, gate_def_targets.args)::Vector{String}
        v.gate_defs[gate_name] = GateDefinition(gate_name, argument_names, qubit_targets, gate_body)
    elseif head(program_expr) == :function_call
        function_name = name(program_expr)
        if haskey(builtin_functions, function_name)
            concrete_arguments = v(convert(Vector{QasmExpression}, program_expr.args[2].args))
            if function_name != "sizeof"
                return_val = builtin_functions[function_name](Iterators.flatten(concrete_arguments)...)
            else
                return_val = builtin_functions[function_name](concrete_arguments...)
            end
            return return_val[1]
        else
            hasfunction(v, function_name) || throw(QasmVisitorError("function $function_name not defined!"))
            function_def  = function_defs(v)[function_name]
            function_body = function_def.body::Vector{QasmExpression}
            declared_args = only(function_def.arguments.args)::QasmExpression
            provided_args = only(program_expr.args[2].args)::QasmExpression
            function_v    = QasmFunctionVisitor(v, declared_args, provided_args)
            return_val    = nothing
            body_exprs::Vector{QasmExpression} = head(function_body[1]) == :scope ? function_body[1].args : function_body
            for f_expr in body_exprs
                if head(f_expr) == :return
                    return_val = function_v(f_expr.args[1])
                else
                    function_v(f_expr)
                end
            end
            # remap qubits and classical variables
            function_args = if head(declared_args) == :array_literal
                convert(Vector{QasmExpression}, declared_args.args)::Vector{QasmExpression}
            else
                declared_args
            end
            called_args = if head(provided_args) == :array_literal
                convert(Vector{QasmExpression}, provided_args.args)::Vector{QasmExpression}
            else
                provided_args
            end
            arguments_map         = Dict{QasmExpression, QasmExpression}(zip(function_args, called_args))
            reverse_arguments_map = Dict{QasmExpression, QasmExpression}(zip(called_args, function_args))
            reverse_qubits_map    = Dict{Int, Int}()
            for variable in keys(reverse_arguments_map)
                if head(variable) ∈ (:identifier, :indexed_identifier)
                    variable_name = name(variable)
                    if haskey(classical_defs(v), variable_name) && classical_defs(v)[variable_name].type isa SizedArray
                        if head(reverse_arguments_map[variable]) != :const_declaration
                            inner_variable_name = name(reverse_arguments_map[variable])
                            new_val = classical_defs(function_v)[inner_variable_name].val
                            back_assignment = QasmExpression(:classical_assignment, QasmExpression(:binary_op, Symbol("="), variable, new_val))
                            v(back_assignment)
                        end
                    elseif haskey(qubit_defs(v), variable_name)
                        outer_context_map = only(evaluate_qubits(v, variable))
                        inner_context_map = only(evaluate_qubits(function_v, reverse_arguments_map[variable].args[1]))
                        reverse_qubits_map[inner_context_map] = outer_context_map
                    end
                end
            end
            mapper = isempty(reverse_qubits_map) ? identity : ix->remap(ix, reverse_qubits_map)
            push!(v, map(mapper, function_v.instructions))
            return return_val
        end
    elseif head(program_expr) == :function_definition
        function_def         = program_expr.args
        function_name        = function_def[1].args[1]::String
        function_arguments   = function_def[2]
        function_return_type = function_def[3]::QasmExpression
        function_body        = function_def[4]::QasmExpression
        full_function_def    = FunctionDefinition(function_name, function_arguments, function_body, function_return_type)
        v.function_defs[function_name] = full_function_def
    elseif head(program_expr) == :pragma
        visit_pragma(v, program_expr)
    elseif head(program_expr) ∈ (:integer_literal, :float_literal, :string_literal, :complex_literal, :irrational_literal, :boolean_literal, :duration_literal)
        return program_expr.args[1]
    elseif head(program_expr) == :array_literal
        return map(v, program_expr.args)
    elseif head(program_expr) == :range
        start::Int, step::Int, stop::Int = v(program_expr.args)
        return StepRange(start, step, stop)
    elseif head(program_expr) == :empty
        return ()
    elseif head(program_expr) == :measure
        qubits_to_measure = evaluate_qubits(v, program_expr.args[1])
        push!(v, CircuitInstruction[(type="measure", arguments=InstructionArgument[], targets=[q], controls=Pair{Int,Int}[], exponent=1.0) for q in qubits_to_measure])
        return false
    elseif head(program_expr) == :hw_qubit
        return tryparse(Int, name(program_expr))
    elseif head(program_expr) == :output
        throw(QasmVisitorError("Output not supported."))
    elseif head(program_expr) == :binary_op
        op  = program_expr.args[1]::Symbol
        lhs = v(program_expr.args[2])
        rhs = v(program_expr.args[3])
        val = evaluate_binary_op(op, lhs, rhs)
        return val
    elseif head(program_expr) == :unary_op
        op  = program_expr.args[1]::Symbol
        arg = v(program_expr.args[2])
        return evaluate_unary_op(op, arg)
    elseif head(program_expr) == :cast
        casting_to = program_expr.args[1].args[1]
        value = v(program_expr.args[2])
        if casting_to == Bool
            return value > 0
        # TODO
        else
            throw(QasmVisitorError("unable to evaluate cast expression $program_expr"))
        end
    else
        throw(QasmVisitorError("cannot visit expression $program_expr.")) 
    end
    return v
end

end # module Quasar
