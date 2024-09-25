function _error_representation(observable::Matrix{ComplexF64})
    mat = Vector{Vector{Vector{Float64}}}(undef, size(observable, 1))
    for row in 1:size(observable, 1)
        mat[row] = Vector{Vector{Float64}}(undef, size(observable, 2))
        for col in 1:size(observable, 2)
            mat[row][col] = [real(observable[row, col]), imag(observable[row, col])]
        end
    end
    return mat
end
_error_representation(observable::String) = observable

function _check_observable_targets(observable::Union{Matrix{ComplexF64}, String}, targets)
    qc = Quasar.qubit_count(observable)
    qc == 1 && (isempty(targets) || length(targets) == 1) && return
    qc == length(targets) && return
    throw(Quasar.QasmVisitorError("Invalid observable specified: $(_error_representation(observable)), targets: $targets", "ValueError"))
end
function _check_observable_targets(observable::Vector{Union{String, Matrix{ComplexF64}}}, targets)
    if length(observable) == 1
        _check_observable_targets(only(observable), targets)
    else
        return nothing
    end
end
_check_observable_targets(observable, targets) = nothing

function visit_observable(v, expr)
    raw_obs = expr.args[1]::Quasar.QasmExpression
    if Quasar.head(raw_obs) == :array_literal
        return mapreduce(arg->visit_observable(v, Quasar.QasmExpression(:observable, arg)), vcat, convert(Vector{Quasar.QasmExpression}, raw_obs.args))
    elseif Quasar.head(raw_obs) == :identifier
        return Union{String, Matrix{ComplexF64}}[raw_obs.args[1]]
    elseif Quasar.head(raw_obs) == :hermitian
        h_mat = similar(raw_obs.args[1], ComplexF64)::Matrix{ComplexF64}
        for ii in eachindex(h_mat)
            h_mat[ii] = convert(ComplexF64, v(raw_obs.args[1][ii]))::ComplexF64
        end
        return Union{String, Matrix{ComplexF64}}[h_mat]
    end
end

function Quasar.visit_pragma(v, program_expr)
    pragma_type::Symbol = program_expr.args[1]
    if pragma_type == :result
        result_type = program_expr.args[2]
        if result_type == :state_vector
            push!(v, (type=:state_vector, operator=Union{String, Matrix{ComplexF64}}[], targets=Int[], states=String[]))
        elseif result_type ∈ (:probability, :density_matrix)
            has_targets = !isempty(program_expr.args[3].args)
            targets = has_targets ? Quasar.evaluate_qubits(v, program_expr.args[3].args[1]) : Int[]
            rt = (type=result_type, operator=Union{String, Matrix{ComplexF64}}[], targets=targets, states=String[])
            push!(v, rt)
        elseif result_type == :amplitude
            states = Quasar.head(program_expr.args[3]) == :array_literal ? program_expr.args[3].args : program_expr.args[3]
            clean_states = map(state->replace(state.args[1], "\""=>"", "\'"=>""), states)
            push!(v, (type=:amplitude, operator=Union{String, Matrix{ComplexF64}}[], targets=Int[], states=clean_states))
        elseif result_type ∈ (:expectation, :variance, :sample)
            raw_obs     = program_expr.args[3]::Quasar.QasmExpression
            raw_targets = program_expr.args[end]::Quasar.QasmExpression
            has_targets = !isempty(raw_targets.args)
            targets     = has_targets ? Quasar.evaluate_qubits(v, raw_targets.args[1]) : Int[] 
            observable  = visit_observable(v, raw_obs)
            if length(observable) == 1 && only(observable) ∈ ("x", "y", "z", "h", "i") && length(targets) > 1
                throw(Quasar.QasmVisitorError("Standard observable target must be exactly 1 qubit.", "ValueError"))
            end
            _check_observable_targets(observable, targets)
            push!(v, (type=result_type, operator=observable, targets=targets, states=String[]))
        elseif result_type == :adjoint_gradient
            throw(Quasar.QasmVisitorError("Result type adjoint_gradient is not supported.", "TypeError"))
        end
    elseif pragma_type == :unitary
        raw_mat = program_expr.args[2]::Matrix{Quasar.QasmExpression}
        unitary_matrix = similar(raw_mat, ComplexF64)::Matrix{ComplexF64}
        for ii in eachindex(unitary_matrix)
            unitary_matrix[ii] = v(raw_mat[ii])
        end
        targets = Quasar.evaluate_qubits(v, program_expr.args[end].args[1])
        push!(v, (type="unitary", arguments=Union{Symbol, Dates.Period, Real, Matrix{ComplexF64}}[unitary_matrix], targets=targets, controls=Pair{Int,Int}[], exponent=1.0))
    elseif pragma_type == :noise
        noise_type::String = program_expr.args[2].args[1]::String
        raw_args::Quasar.QasmExpression = program_expr.args[3].args[1]::Quasar.QasmExpression
        raw_targets::Quasar.QasmExpression = program_expr.args[4]::Quasar.QasmExpression
        targets = Quasar.evaluate_qubits(v, raw_targets.args[1])::Vector{Int}
        if noise_type == "kraus"
            raw_mats = raw_args.args
            kraus_matrices = Union{Symbol, Dates.Period, Real, Matrix{ComplexF64}}[broadcast(expr->convert(ComplexF64, v(expr)), raw_mat)::Matrix{ComplexF64} for raw_mat in raw_mats]
            push!(v, (type="kraus", arguments=kraus_matrices, targets=targets, controls=Pair{Int,Int}[], exponent=1.0))
        else
            args = map(Float64, v(raw_args))
            push!(v, (type=noise_type, arguments=Union{Symbol, Dates.Period, Real, Matrix{ComplexF64}}[args...], targets=targets, controls=Pair{Int,Int}[], exponent=1.0))
        end
    elseif pragma_type == :verbatim
    end
    return v
end

function parse_pragma_observables(tokens::Vector{Tuple{Int64, Int32, Quasar.Token}}, stack, start, qasm)
    observables_list = Quasar.QasmExpression[]
    obs_targets      = Quasar.QasmExpression[]
    while !isempty(tokens)
        observable_token = popfirst!(tokens)
        observable_id    = Quasar.parse_identifier(observable_token, qasm)
        if observable_id.args[1] == "hermitian"
            matrix_tokens = Quasar.extract_parensed(tokens, stack, start, qasm)
            # next token is targets
            h_mat = Quasar.parse_matrix(matrix_tokens, stack, start, qasm)
            # next token is targets
            next_at = findfirst(triplet->triplet[end] == Quasar.at, tokens)
            final_token = isnothing(next_at) ? length(tokens) : next_at-1
            target_tokens = splice!(tokens, 1:final_token)
            if !(isempty(target_tokens) || first(target_tokens)[end] == Quasar.all_token)
                push!(target_tokens, (-1, Int32(-1), Quasar.semicolon))
                while !isempty(target_tokens) && first(target_tokens)[end] != Quasar.semicolon
                    target_tokens[1][end] == Quasar.comma && popfirst!(target_tokens)
                    target_expr = Quasar.parse_expression(target_tokens, stack, target_tokens[1][1], qasm)
                    push!(obs_targets, target_expr)
                end
            end
            push!(observables_list, Quasar.QasmExpression(:hermitian, h_mat))
        elseif observable_id.args[1] == "all"
            break
        else
            if !isempty(tokens) && first(tokens)[end] == Quasar.lparen
                arg_tokens = Quasar.extract_parensed(tokens, stack, start, qasm)
                push!(arg_tokens, (-1, Int32(-1), Quasar.semicolon))
                target_expr = Quasar.parse_expression(arg_tokens, stack, start, qasm)
                push!(obs_targets, target_expr)
            end
            push!(observables_list, observable_id)
        end
        !isempty(tokens) && first(tokens)[end] == Quasar.at && popfirst!(tokens)
    end
    if length(observables_list) == 1 && length(obs_targets) == 1
        return Quasar.QasmExpression(:observable, only(observables_list)), Quasar.QasmExpression(:qubit_targets, only(obs_targets))
    elseif length(observables_list) == 1 && length(obs_targets) == 0
        return Quasar.QasmExpression(:observable, only(observables_list)), Quasar.QasmExpression(:qubit_targets)
    elseif length(observables_list) == 1 && length(obs_targets) > 1
        return Quasar.QasmExpression(:observable, only(observables_list)), Quasar.QasmExpression(:qubit_targets, Quasar.QasmExpression(:array_literal, obs_targets...))
    else
        return Quasar.QasmExpression(:observable, Quasar.QasmExpression(:array_literal, observables_list...)), Quasar.QasmExpression(:qubit_targets, Quasar.QasmExpression(:array_literal, obs_targets...))
    end
end

function parse_pragma_targets(tokens::Vector{Tuple{Int64, Int32, Quasar.Token}}, stack, start, qasm)
    target_expr = Quasar.QasmExpression(:qubit_targets)
    (isempty(tokens) || first(tokens)[end] == Quasar.all_token) && return target_expr
    push!(tokens, (-1, Int32(-1), Quasar.semicolon))
    push!(target_expr, Quasar.parse_list_expression(tokens, stack, start, qasm))
    return target_expr
end


function Quasar.parse_pragma(tokens, stack, start, qasm)
    prefix    = popfirst!(tokens)
    prefix_id = Quasar.parse_identifier(prefix, qasm)
    prefix_id.args[1] == "braket" || throw(Quasar.QasmParseError("pragma expression must begin with `#pragma braket`", stack, start, qasm))
    expr        = Quasar.QasmExpression(:pragma)
    pragma_type = Quasar.parse_identifier(popfirst!(tokens), qasm).args[1]
    if pragma_type == "result"
        push!(expr, :result)
        result_type = Quasar.parse_identifier(popfirst!(tokens), qasm).args[1]
        if result_type == "state_vector"
            push!(expr, :state_vector)
        elseif result_type == "probability"
            target_expr = parse_pragma_targets(tokens, stack, start, qasm)
            push!(expr, :probability, target_expr)
        elseif result_type == "density_matrix"
            target_expr = parse_pragma_targets(tokens, stack, start, qasm)
            push!(expr, :density_matrix, target_expr)
        elseif result_type == "amplitude"
            push!(tokens, (-1, Int32(-1), Quasar.semicolon))
            states = Quasar.parse_list_expression(tokens, stack, start, qasm)
            push!(expr, :amplitude, states)
        elseif result_type ∈ ("expectation", "variance", "sample")
            obs, targets = parse_pragma_observables(tokens, stack, start, qasm)
            push!(expr, Symbol(result_type), obs, targets)
        elseif result_type == "adjoint_gradient"
            push!(expr, :adjoint_gradient)
        end
    elseif pragma_type == "unitary"
        push!(expr, :unitary)
        matrix_tokens  = Quasar.extract_parensed(tokens, stack, start, qasm)
        unitary_matrix = Quasar.parse_matrix(matrix_tokens, stack, start, qasm)
        push!(expr, unitary_matrix)
        target_expr = parse_pragma_targets(tokens, stack, start, qasm)
        push!(expr, target_expr)
    elseif pragma_type == "noise"
        push!(expr, :noise)
        noise_type = Quasar.parse_identifier(popfirst!(tokens), qasm)::Quasar.QasmExpression
        if noise_type.args[1] == "kraus"
            matrix_tokens = Quasar.extract_parensed(tokens, stack, start, qasm)
            all(triplet->triplet[end] == Quasar.lbracket, matrix_tokens[1:3]) && (matrix_tokens = Quasar.extract_braced_block(matrix_tokens, stack, start, qasm))
            mats = Matrix{Quasar.QasmExpression}[]
            while !isempty(matrix_tokens)
                push!(mats, Quasar.parse_matrix(matrix_tokens, stack, start, qasm))
                isempty(matrix_tokens) && break
                next_token = first(matrix_tokens)
                next_token[end] == Quasar.comma && popfirst!(matrix_tokens)
                next_token[end] == Quasar.semicolon && break
            end
            noise_args = Quasar.QasmExpression(:arguments, Quasar.QasmExpression(:array_literal, mats))
        else
            noise_args = Quasar.parse_arguments_list(tokens, stack, start, qasm)
        end
        push!(expr, noise_type, noise_args)
        target_expr = parse_pragma_targets(tokens, stack, start, qasm)
        push!(expr, target_expr)
    elseif pragma_type == "verbatim"
        # check that the next non-newline is a box token 
        push!(expr, :verbatim)
    else
        throw(Quasar.QasmParseError("invalid type $pragma_type for pragma", stack, start, qasm))
    end
    return expr
end
