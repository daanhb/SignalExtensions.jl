# sequences.jl

"""
A Sequence is a bi-infinite discrete set of values, indexed by all integers.

Each Sequence has an element type, given by `eltype`.
"""
abstract Sequence{T}

Base.eltype{T}(::Type{Sequence{T}}) = T
Base.eltype{S <: Sequence}(::Type{S}) = eltype(supertype(S))

# Return a range of the sequence as a vector
Base.getindex(s::Sequence, r::Range) = [s[k] for k in range]

iscompact(s::Sequence) = false

# We adopt the convention that the transpose of a sequence corresponds to time-reversal
Base.transpose(s::Sequence) = reverse(s)

# Hence, ctranspose is time-reversal plus taking conjugates
Base.ctranspose(s::Sequence) = conj(reverse(s))
