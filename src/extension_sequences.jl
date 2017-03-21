# extensionsequences.jl


#################################################
# Abstract extension type and generic interface
#################################################

# In general, extension sequences specialize on the type of the underlying
# data. This means operations with extension sequences are typically fast.
# Since a point to an array is stored, memory is allocated every time an
# extension sequence is constructed. This may be avoided in future versions
# of Julia.

"""
Any subtype of ExtensionSequence embeds an indexable vectorlike object with finite length
and extends it into a bi-infinite sequence. The extension can be indexed with any
integer.

Indexing with an index in the allowed range of the embedded vector simply returns
the corresponding element. Indexing with indices outside this range result in
a computation. For example, for periodic extensions, the indices are mapped to
the allowed range using periodicity.

Usually the indices of the embedded vector `a` range from `1` to `length(a)`, but
that need not be the case. The extension respects the range returned by
`linearindices(a)`. One could use an OffsetArray with indices `-2:2`.

In current versions of Julia, construction of an ExtensionSequence allocates
memory, because the types store a reference to the embedded vector.
"""
abstract ExtensionSequence{T} <: Sequence{T}


# We assume that the embedded vector is in the field 'a'.
# A convention we follow is that variable 'k' refers to indices of the sequence,
# and variable 'i' to indices of the embedded vector.

"The embedded vector of the extension sequence."
subvector(s::ExtensionSequence) = s.a

"The length of the embedded vector of the extension sequence."
sublength(s::ExtensionSequence) = length(linearindices(s.a))

"The first linear index of the embedded vector."
first_subindex(s::ExtensionSequence) = first(linearindices(s.a))

"The last linear index of the embedded vector."
last_subindex(s::ExtensionSequence) = last(linearindices(s.a))

# The indexing logic can refer to the `element` function defined for all extension types
# Concrete subtypes of ExtensionSequence have to implement `extensiontype`.
getindex(s::ExtensionSequence, k) = element(extensiontype(s), subvector(s), k)

# By implementing `conj` and `reverse`, we get `transpose` and `ctranspose` for free
# Concrete types should implement a `similar` function that returns a sequence of
# the same type but with a different embedded vector.
for op in (:conj,)
    @eval $op(s::ExtensionSequence) = similar(s, $op(subvector(s)))
end

reverse(s::ExtensionSequence) = _reverse(s, subvector(s))
_reverse(s, a) = similar(s, reverse_indices(a))


#######################
# Periodic extensions
#######################

"""
A PeriodicSequence extends a vector `a` periodically to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range are mapped
by periodization modulo the length of `a`.
"""
immutable PeriodicSequence{A,T} <: ExtensionSequence{T}
    a       :: A
end

PeriodicSequence(a::AbstractVector) = PeriodicSequence{typeof(a),eltype(a)}(a)

similar(p::PeriodicSequence, a) = PeriodicSequence(a)

extensiontype(s::PeriodicSequence) = PeriodicExtension()



#######################
# CompactSequence
#######################

"""
CompactSequence extends a vector `a` with zeros to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range correspond
to zero values.
"""
immutable CompactSequence{A,T} <: ExtensionSequence{T}
    a :: A
end

typealias ZeroPaddingSequence CompactSequence

CompactSequence{A}(a::A) = CompactSequence{A,eltype(A)}(a)

similar(s::CompactSequence, a) = CompactSequence(a)

extensiontype(s::CompactSequence) = ZeroPadding()

iscompact(::CompactSequence) = true

nzrange(s::CompactSequence) = linearindices(subvector(s))


#######################
# Constant padding
#######################


"""
A ConstantPaddingSequence extends a vector `a` with a given `constant` to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range correspond to the constant
value.
"""
immutable ConstantPaddingSequence{A,T} <: ExtensionSequence{T}
    a           ::  A
    constant    ::  T
end

ConstantPaddingSequence(a::AbstractVector, constant) = ConstantPaddingSequence{typeof(a),eltype(a)}(a, constant)

constant(s::ConstantPaddingSequence) = s.constant

similar(s::ConstantPaddingSequence, a) = ConstantPaddingSequence(a, constant(s))

extensiontype(s::ConstantPaddingSequence) = ConstantPadding(s.constant)



#######################
# Symmetric extensions
#######################


"""
A SymmetricSequence extends a vector 'a' symmetrically to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range are mapped by symmetrizing.

The symmetry around each of the endpoints can be whole-point (the endpoint is not repeated)
or half-point (the endpoint is repeated). The symmetry can also be even (symmetric)
or odd (anti-symmetric).

Definition:

immutable SymmetricSequence{SYM_LEFT,SYM_RIGHT,PAR_LEFT,PAR_RIGHT,A,T} <: ExtensionSequence{T}

Type parameters indicate the type of symmetry and its parity:
- SYM_LEFT:  either :wp (whole point) or :hp (half point) near left endpoint
- SYM_RIGHT: either :wp or :hp for right endpoint
- PAR_LEFT: either :odd or :even symmetry near left endpoint
- PAR_RIGHT: also :odd or :even
- A: the type of the embedded vector
- T: the eltype of the sequence
"""
immutable SymmetricSequence{SYM_LEFT,SYM_RIGHT,PAR_LEFT,PAR_RIGHT,A,T} <: ExtensionSequence{T}
    a :: A
end

SymmetricSequence(a::AbstractVector) = symmetric_extension_wholepoint_even(a)

similar{SL,SR,PL,PR}(s::SymmetricSequence{SL,SR,PL,PR}, a) =
    SymmetricSequence{SL,SR,PL,PR,typeof(a),eltype(a)}(a)

# Provide four of the sixteen combinations for convenience. Construct the others by explicitly
# specifying the type parameters symbols.
symmetric_extension_wholepoint_even(a) =
    SymmetricSequence{:wp,:wp,:even,:even,typeof(a),eltype(a)}(a)

symmetric_extension_halfpoint_even(a) =
    SymmetricSequence{:hp,:hp,:even,:even,typeof(a),eltype(a)}(a)

symmetric_extension_wholepoint_odd(a) =
    SymmetricSequence{:wp,:wp,:odd,:odd,typeof(a),eltype(a)}(a)

symmetric_extension_halfpoint_odd(a) =
    SymmetricSequence{:hp,:hp,:odd,:odd,typeof(a),eltype(a)}(a)

# We have to be careful when reversing a symmetric extension: the endpoints change places
function reverse{SL,SR,PL,PR}(s::SymmetricSequence{SL,SR,PL,PR})
    b = reverse_indices(subvector(s))
    SymmetricSequence{SR,SL,PR,PL,typeof(b),eltype(b)}(b)
end

for op in (:left_parity, :right_parity, :left_symmetry, :right_symmetry)
    @eval $op(s::SymmetricSequence) = $op(extensiontype(s))
end

extensiontype{SL,SR,PL,PR}(s::SymmetricSequence{SL,SR,PL,PR}) =
    SymmetricExtension{SL,SR,PL,PR}()
