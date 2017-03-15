# sequences.jl

"""
A Sequence is a bi-infinite discrete set of values, indexed by all integers.

Each Sequence has an element type, given by `eltype`.
"""
abstract Sequence{T}

Base.eltype{T}(::Type{Sequence{T}}) = T

"Return the index of the first non-zero element in the sequence."
firstindex(s::Sequence) = -inf

"Return the index of the last non-zero element in the sequence."
lastindex(s::Sequence) = inf

# By convention, eachindex sums only over nonzero elements of the sequence.
eachindex(s::Sequence) = firstindex(s):lastindex(s)

# Range of indices that covers all nonzero elements of both sequences.
eachindex(s1::Sequence, s2::Sequence) =
    min(firstindex(s1),firstindex(s2)):max(lastindex(s1),lastindex(s2))

"An iterator over all non-zero values in the sequence."
each_nonzero_index(s::Sequence) = eachindex(s)

# Return a range of the sequence as a vector
getindex(s::Sequence, r::Range) = [s[k] for k in range]

hascompactsupport(s::Sequence) = false

transpose(s::Sequence) = reverse(s)

"The j-th discrete moment of a sequence is defined as `\sum_k h_k k^j`."
function moment(s::Sequence, j)
    z = zero(eltype(s))
    for k in each_nonzero_index(s)
        z += s[k] * k^j
    end
    z
end
