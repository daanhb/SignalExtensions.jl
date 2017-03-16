# extensionsequences.jl


#################################################
# Abstract extension type and generic interface
#################################################

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

The ExtensionSequence acts as a mutating view. One can set elements of the
extension, and the corresponding entry of the subvector will be modified.
"""
abstract ExtensionSequence{T} <: Sequence{T}


# We assume that the embedded vector is in the field 'a'.
# A convention we follow is that variable 'k' refers to indices of the sequence,
# and variable 'i' to indices of the embedded vector.
#
# All subtypes should implement:
# mapindex(s::SomeSubType, k) -> map an index 'k' of the sequence to an index 'i'
# of the embedded vector
# imapindex(s::SomeSubType, i) -> the inverse map

"The embedded vector of the extension sequence."
subvector(s::ExtensionSequence) = s.a

"The length of the embedded vector of the extension sequence."
sublength(s::ExtensionSequence) = length(linearindices(s.a))

"The first linear index of the embedded vector."
first_subindex(s::ExtensionSequence) = first(linearindices(s.a))

"The last linear index of the embedded vector."
last_subindex(s::ExtensionSequence) = last(linearindices(s.a))


"Map an index into the extension sequence to an index into the embedded vector."
mapindex(s::ExtensionSequence, k) = k

"Map an index into the embedded vector into an index into the extension sequence."
imapindex(s::ExtensionSequence, i) = i

getindex(s::ExtensionSequence, k) = getindex(s.a, mapindex(s, k))

setindex!(s::ExtensionSequence, val, k) = setindex!(s.a, val, mapindex(s, k))

for op in (:conj,)
    @eval $op(s::ExtensionSequence) = similar(s, $op(subvector(s)))
end

reverse(s::ExtensionSequence) = _reverse(s, subvector(s))
_reverse(s, a) = similar(s, reverse_indices(a))


#######################
# Periodic extensions
#######################

"""
A PeriodicExtension extends a vector `a` periodically to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range are mapped
by periodization modulo the length of `a`.

The periodic extension acts as a mutating view. One can set elements of the
extension, and the corresponding entry of `a` will be modified.
"""
immutable PeriodicExtension{A,T} <: ExtensionSequence{T}
    a       :: A
end

PeriodicExtension(a::AbstractVector) = PeriodicExtension{typeof(a),eltype(a)}(a)

similar(p::PeriodicExtension, a) = PeriodicExtension(a)

mapindex(s::PeriodicExtension, k) = _mapindex(s, k, first_subindex(s), last_subindex(s), sublength(s))

# Note that we only do the expensive modulo operator after the bounds checking
_mapindex(s::PeriodicExtension, k, i0, i1, n) = i0 <= k <= i1 ? k : mod(k-i0, n) + i0

# No need to override imapindex


#######################
# ZeroPadding
#######################

"""
ZeroPadding extends a vector `a` with zeros to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range correspond
to zero values.

A ZeroPadding acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of `a` will be modified.
"""
immutable ZeroPadding{A,T} <: ExtensionSequence{T}
    a :: A
end

ZeroPadding{A}(a::A) = ZeroPadding{A,eltype(A)}(a)

similar(s::ZeroPadding, a) = ZeroPadding(a)

# We override getindex to return zero outside of our embedded vector.
getindex(s::ZeroPadding, k::Int) = (k < first_subindex(s)) || (k > last_subindex(s)) ? zero(eltype(s)) : getindex(s.a, k)

# No need to override mapindex and imapindex

iscompact(::ZeroPadding) = true


#######################
# Constant padding
#######################


"""
A ConstantPadding extends a vector `a` with a given `constant` to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range correspond to the constant
value.

A ConstantPadding acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of `a` will be modified.
"""
immutable ConstantPadding{A,T} <: ExtensionSequence{T}
    a           ::  A
    constant    ::  T
end

ConstantPadding(a::AbstractVector, constant) = ConstantPadding{typeof(a),eltype(a)}(a, constant)

constant(s::ConstantPadding) = s.constant

similar(s::ConstantPadding, a) = ConstantPadding(a, constant(s))

# We override getindex to return the constant outside of our embedded vector.
getindex(s::ConstantPadding, k::Int) = (k < first_subindex(s)) || (k > last_subindex(s)) ? s.constant : getindex(s.a, k)

# No need to override mapindex and imapindex


#######################
# Symmetric extensions
#######################


"""
A SymmetricExtension extends a vector 'a' symmetrically to a bi-infinite sequence.

The indices of `a` map to `a`. Indices outside this range are mapped by symmetrizing.

The symmetry around each of the endpoints can be whole-point (the endpoint is not repeated)
or half-point (the endpoint is repeated). The symmetry can also be even (symmetric)
or odd (anti-symmetric).

The symmetric extension acts as a mutating view. One can set elements of the
extension, and the corresponding entry of 'a' will be modified.

Definition:

immutable SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT} <: ExtensionSequence{A}

Type parameters:
- A:    the type of the embedded vector
- PT_LEFT:  either :wp (whole point) or :hp (half point) near left endpoint
- PT_RIGHT: either :wp or :hp for right endpoint
- SYM_LEFT: either :odd or :even symmetry near left endpoint
- SYM_RIGHT: also :odd or :even

"""
immutable SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT,A,T} <: ExtensionSequence{T}
    a :: A
end

SymmetricExtension(a::AbstractVector) = symmetric_extension_wholepoint_even(a)

similar{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}, a) =
    SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT,typeof(a),eltype(a)}(a)

# Provide four of the sixteen combinations for convenience. Construct the others by explicitly
# specifying the type parameters symbols.
symmetric_extension_wholepoint_even(a) =
    SymmetricExtension{:wp,:wp,:even,:even,typeof(a),eltype(a)}(a)

symmetric_extension_halfpoint_even(a) =
    SymmetricExtension{:hp,:hp,:even,:even,typeof(a),eltype(a)}(a)

symmetric_extension_wholepoint_odd(a) =
    SymmetricExtension{:wp,:wp,:odd,:odd,typeof(a),eltype(a)}(a)

symmetric_extension_halfpoint_odd(a) =
    SymmetricExtension{:hp,:hp,:odd,:odd,typeof(a),eltype(a)}(a)

# We have to be careful when reversing a symmetric extension: the endpoints change places
function reverse{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})
    b = reverse_indices(subvector(s))
    SymmetricExtension{PT_RIGHT,PT_LEFT,SYM_RIGHT,SYM_LEFT,typeof(b),eltype(b)}(b)
end

left_parity{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})    = SYM_LEFT
right_parity{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})   = SYM_RIGHT
left_symmetry{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})  = PT_LEFT
right_symmetry{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}) = PT_RIGHT


# Compute the index by mapping any index outside the range of the embedded vector to
# an index that is closer to the interval and repeat. The recursion ends when the index
# lands inside the interval.
# For far away points this is not the fastest way, but that won't happen very often. The
# alternative is to reduce the extension by exploiting periodicity of the symmetric extension,
# but the period depends on all symmetries too and this will require a modulo operation. For
# nearby points, the implementation below is probably faster.
function mapindex(s::SymmetricExtension, k)
    if k > last_subindex(s)
        # We are on the right of the interval: use symmetry wrt right endpoint
        mapindex_right(s, k)
    elseif k < first_subindex(s)
        # We are on the left of the interval: use symmetry wrt left endpoint
        mapindex_left(s, k)
    else
        k
    end
end

# Right whole point symmetry:
mapindex_right{PT_LEFT}(s::SymmetricExtension{PT_LEFT,:wp}, k) = mapindex(s, 2*last_subindex(s) - k + 1)

# Right half point symmetry:
mapindex_right{PT_LEFT}(s::SymmetricExtension{PT_LEFT,:hp}, k) = mapindex(s, 2*last_subindex(s) - k)

# Left whole point symmetry:
mapindex_left{PT_RIGHT}(s::SymmetricExtension{:wp,PT_RIGHT}, k) = mapindex(s, 2*first_subindex(s) - k - 1)

# Left half point symmetry:
mapindex_left{PT_RIGHT}(s::SymmetricExtension{:hp,PT_RIGHT}, k) = mapindex(s, 2*first_subindex(s) - k)

# No need to override imapindex

# For getindex we have to use the same logic as mapindex, but now we also have to trace
# the odd/even-ness of the symmetries.
function getindex(s::SymmetricExtension, k)
    if k > last_subindex(s)
        # We are on the right of the interval: use symmetry wrt right endpoint
        getindex_right(s, k)
    elseif k < first_subindex(s)
        # We are on the left of the interval: use symmetry wrt left endpoint
        getindex_left(s, k)
    else
        getindex(s.a, k)
    end
end

# For even symmetry on both endpoints, it is simple: the sign never flips. Short circuit.
# This probably does not gain much, as mapindex still does the recursion anyway...
getindex{PT_LEFT,PT_RIGHT}(s::SymmetricExtension{PT_LEFT,PT_RIGHT,:even,:even}, k) = getindex(s.a, mapindex(s, k))

# Right whole point symmetry
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricExtension{PT_LEFT,:wp,SYM_LEFT,:even}, k) = getindex(s, 2*last_subindex(s) - k + 1)
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricExtension{PT_LEFT,:wp,SYM_LEFT,:odd}, k) = -getindex(s, 2*last_subindex(s) - k + 1)

# Right half point symmetry
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricExtension{PT_LEFT,:hp,SYM_LEFT,:even}, k) = getindex(s, 2*last_subindex(s) - k)
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricExtension{PT_LEFT,:hp,SYM_LEFT,:odd}, k) = -getindex(s, 2*last_subindex(s) - k)

# Left whole point symmetry
getindex_left{PT_RIGHT}(s::SymmetricExtension{:wp,PT_RIGHT,:even}, k) = getindex(s, 2*first_subindex(s) - k - 1)
getindex_left{PT_RIGHT}(s::SymmetricExtension{:wp,PT_RIGHT,:odd}, k) = -getindex(s, 2*first_subindex(s) - k - 1)

# Left half point symmetry
getindex_left{PT_RIGHT}(s::SymmetricExtension{:hp,PT_RIGHT,:even}, k) = getindex(s, 2*first_subindex(s) - k)
getindex_left{PT_RIGHT}(s::SymmetricExtension{:hp,PT_RIGHT,:odd}, k) = -getindex(s, 2*first_subindex(s) - k)
