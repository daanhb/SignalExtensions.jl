# extensions.jl

abstract Extension

immutable PeriodicExtension <: Extension
end

immutable ZeroPadding <: Extension
end

immutable ConstantPadding{T} <: Extension
    constant    ::  T
end

immutable SymmetricExtension{SYM_LEFT,SYM_RIGHT,PAR_LEFT,PAR_RIGHT} <: Extension
end

left_parity{SL,SR,PL,PR}(s::SymmetricExtension{SL,SR,PL,PR})    = PL
right_parity{SL,SR,PL,PR}(s::SymmetricExtension{SL,SR,PL,PR})   = PR
left_symmetry{SL,SR,PL,PR}(s::SymmetricExtension{SL,SR,PL,PR})  = SL
right_symmetry{SL,SR,PL,PR}(s::SymmetricExtension{SL,SR,PL,PR}) = SR



first_index(a) = first(linearindices(a))
last_index(a) = last(linearindices(a))

"""
Return element `k` from the extension sequence induced by `a`.
"""
element(s::Extension, a, k) = _element(s, a, k, first_index(a), last_index(a))

# Note that we do a boundscheck for `k` here. If `k` is a valid index for `a` we
# use it right away, avoiding other (more expensive) index calculations.
_element(s::Extension, a, k, i0, i1) = (i0 <= k <= i1) ? a[k] : element_extension(s, a, k, i0, i1)


# For periodic functions, extension amounts to a modulo operation.
periodic_index(k, i0, i1) = mod(k-i0, i1-i0+1)+i0

element_extension(s::PeriodicExtension, a, k, i0, i1) = a[periodic_index(k, i0, i1)]

# For zero padding of course we return zero (of the same type as elements of a)
element_extension(s::ZeroPadding, a, k, i0, i1) = zero(eltype(a))

# For constant padding we return the stored constant, converted to eltype(a)
element_extension(s::ConstantPadding, a, k, i0, i1) = convert(eltype(a), s.constant)


# The symmetric case is a little more involved.
#
# Methodology: wenever the index `k` is outside the allowed range, we fold it back
# using symmetry and we iterate. With each iteration, it goes closer to the allowed range.
# An alternative for large indices could be based on the periodicity of the symmetric
# extension.

function _element(s::SymmetricExtension, a, k, i0, i1)
    if k > i1
        # We are to the right of the range i0:i1
        element_right(s, a, k, i0, i1)
    elseif k < i0
        # We are to the left of the range i0:i1
        element_left(s, a, k, i0, i1)
    else
        # k is a valid index for a
        a[k]
    end
end

# Right whole point symmetry, even or odd
element_right{SYM_LEFT,PAR_LEFT}(s::SymmetricExtension{SYM_LEFT,:wp,PAR_LEFT,:even}, a, k, i0, i1) = _element(s, a, 2*i1 - k + 1, i0, i1)
element_right{SYM_LEFT,PAR_LEFT}(s::SymmetricExtension{SYM_LEFT,:wp,PAR_LEFT,:odd}, a, k, i0, i1) = -_element(s, a, 2*i1 - k + 1, i0, i1)

# Right half point symmetry, even or odd
element_right{SYM_LEFT,PAR_LEFT}(s::SymmetricExtension{SYM_LEFT,:hp,PAR_LEFT,:even}, a, k, i0, i1) = _element(s, a, 2*i1 - k, i0, i1)
element_right{SYM_LEFT,PAR_LEFT}(s::SymmetricExtension{SYM_LEFT,:hp,PAR_LEFT,:odd}, a, k, i0, i1) = -_element(s, a, 2*i1 - k, i0, i1)

# Left whole point symmetry, even or odd
element_left{SYM_RIGHT}(s::SymmetricExtension{:wp,SYM_RIGHT,:even}, a, k, i0, i1) = _element(s, a, 2*i0 - k - 1, i0, i1)
element_left{SYM_RIGHT}(s::SymmetricExtension{:wp,SYM_RIGHT,:odd}, a, k, i0, i1) = -_element(s, a, 2*i0 - k - 1, i0, i1)

# Left half point symmetry, even or odd
element_left{SYM_RIGHT}(s::SymmetricExtension{:hp,SYM_RIGHT,:even}, a, k, i0, i1) = _element(s, a, 2*i0 - k, i0, i1)
element_left{SYM_RIGHT}(s::SymmetricExtension{:hp,SYM_RIGHT,:odd}, a, k, i0, i1) = -_element(s, a, 2*i0 - k, i0, i1)
