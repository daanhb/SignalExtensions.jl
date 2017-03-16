# sequences.jl

"""
A Sequence is a bi-infinite discrete set of values, indexed by all integers.

Each Sequence has an element type, given by `eltype`.
"""
abstract Sequence{T}

Base.eltype{T}(::Type{Sequence{T}}) = T

# Return a range of the sequence as a vector
Base.getindex(s::Sequence, r::Range) = [s[k] for k in range]

iscompact(s::Sequence) = false

Base.transpose(s::Sequence) = reverse(s)
