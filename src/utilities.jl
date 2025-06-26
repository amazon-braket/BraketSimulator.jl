"""
	pad_bit(i::Int, pos::Int) -> Int

Insert a 0 bit at position pos in the binary representation of i.
This is used for bit manipulation when applying gates to qubits.
"""

function pad_bit(i::Int, pos::Int)
	# Create a mask for bits before the position
	mask_before = (1 << pos) - 1
	# Create a mask for bits after the position
	mask_after = ~mask_before

	# Extract bits before and after the position
	bits_before = i & mask_before
	bits_after = i & mask_after

	# Shift bits after the position to make room for the new bit
	bits_after = bits_after << 1

	# Combine the bits
	return bits_before | bits_after
end

"""
	flip_bit(i::Int, pos::Int) -> Int

Flip the bit at position pos in the binary representation of i.
This is used for bit manipulation when applying gates to qubits.
"""
function flip_bit(i::Int, pos::Int)
	# Create a mask with a 1 at the position to flip
	mask = 1 << pos

	# XOR with the mask to flip the bit
	return i ⊻ mask
end

"""
	evaluate_binary_op(op::Symbol, lhs, rhs)

Evaluate binary operations. This is a direct copy of the visitor's function.
"""

function evaluate_binary_op(op::Symbol, lhs, rhs)
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
	error("Unknown binary operator: op")
end

"""
	evaluate_unary_op(op::Symbol, arg)

Evaluate unary operations. This is a direct copy of the visitor's function.
"""
function evaluate_unary_op(op::Symbol, arg)
	op == :! && return !arg
	op == :~ && return .!arg
	op == :- && return -arg
	error("Unknown unary operator: op")
end