# derived_sequence.jl

"""
A DerivedSequence inherits functionality from an underlying sequence that is
stored. Yet, it is its own type, and hence specific functionality can be added.

Any type that inherits from DerivedSequence has a full interface of a Sequence.
Without any other definitions, the new type is functionally equivalent to the
underlying one.
"""
abstract DerivedSequence{T} <: Sequence{T}
end

"Return the underlying sequence of the derived sequence."
supersequence(s::DerivedSequence) = s.supersequence
# We assume that the underlying sequence is stored in a field called supersequence.

for op in (:iscompact, :nzrange, )
    @eval $op(s::DerivedSequence) = $op(supersequence(s))
end


"An example of a concrete derived sequence."
immutable ConcreteDerivedSequence{T} <: Sequence{T}
    supersequence   ::  Sequence{T}
end
