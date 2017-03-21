# lazy.jl

# A collection of lazy sequences that perform a computation on-the-fly whenever
# getindex is called.
#
# In general, these lazy sequences do not specialize on the type of the underlying
# sequence(s). This avoids nested type parameters, yet it means operations
# are typically slow.

# A lazy sequence can be transformed into a real one using `collect`
function collect(s::Sequence)
    @assert iscompact(s)
    r = nzrange(s)
    data = OffsetArray([s[i] for i in r], r)
    ZeroPaddingSequence(data)
end


"""
A lazy representation of the convolution of two sequences `f` and `g`. Their
convolution `h = f * g` is defined element-wise by:
`h[n] = \sum_k f[k]*g[n-k].`

The lazy convolution performs the above calculation on the fly with every
call to `getindex`.
"""
immutable Convolution{T} <: Sequence{T}
    f   ::  Sequence{T}
    g   ::  Sequence{T}
end

# * for sequences means convolution
(*)(f::Sequence, g::Sequence) = convolve(f, g)

convolve(f::Sequence, g::Sequence) = Convolution(f, g)

getindex(s::Convolution, n) = _getindex(s, n, s.f, s.g)

function _getindex(s::Convolution, n, f, g)
    z = zero(eltype(s))
    if iscompact(f)
        for k in nzrange(f)
            z += f[k] * g[n-k]
        end
    elseif iscompact(g)
        for k in nzrange(g)
            z += g[k] * f[n-k]
        end
    else
        throw(InexactError())
    end
    z
end

iscompact(s::Convolution) = _iscompact(s, s.f, s.g)
_iscompact(s::Convolution, f, g) = iscompact(f) && iscompact(g)

nzrange(s::Convolution) = _nzrange(s, s.f, s.g)

function _nzrange(s::Convolution, f, g)
    @assert _iscompact(s, f, g)
    # We have to determine the range of indices such that f[k] and g[n-k] overlap.
    a = first(nzrange(f))
    b = last(nzrange(f))
    c = first(nzrange(g))
    d = last(nzrange(g))
    # The intervals are [a,b] and [n-d,n-c].
    # Hence, the lower bound is when n-c = a.
    # The upper bound is when n-d = b
    a+c:d+b
end

ztransform(s::Convolution, z) = _ztransform(s, z, s.f, s.g)
_ztransform(s::Convolution, z, f, g) = ztransform(f, z) * ztransform(g, z)


"""
A ShiftedSequence is the lazy representation of a sequence that is shifted by
`k` units:
`g[n] = x[n-k]`.
"""
immutable ShiftedSequence{T} <: Sequence{T}
    f   ::  Sequence{T}
    k   ::  Int
end

sequence_shift(s::ShiftedSequence) = s.k

"Shift a sequence by `k` positions forward."
shift(s::Sequence, k) = ShiftedSequence(s, k)

# Don't shift a shifted sequence twice
shift(s::ShiftedSequence, k) = ShiftedSequence(s.f, s.k+k)

getindex(s::ShiftedSequence, n) = _getindex(s, n, s.f, s.k)
_getindex(s::ShiftedSequence, n, f, k) = f[n-k]

iscompact(s::ShiftedSequence) = _iscompact(s, s.f)
_iscompact(s::ShiftedSequence, f) = iscompact(f)

nzrange(s::ShiftedSequence) = _nzrange(s, s.f, s.k)

function _nzrange(s::ShiftedSequence, f, k)
    @assert _iscompact(s, f)
    r = nzrange(f)
    first(r)+k:last(r)+k
end

ztransform(s::ShiftedSequence, z) = _ztransform(s, z, s.f, s.k)
_ztransform(s::ShiftedSequence, z, f, k) = z^(-k)*ztransform(f, z)


"""
A ReversedSequence is the lazy representation of a sequence that is reversed in
time:
`g[k] = x[-k]`.
"""
immutable ReversedSequence{T} <: Sequence{T}
    f   ::  Sequence{T}
end

reverse(s::Sequence) = ReversedSequence(s)
reverse(s::ReversedSequence) = s.f

getindex(s::ReversedSequence, k) = _getindex(s, k, s.f)
_getindex(s::ReversedSequence, k, f) = f[-k]

iscompact(s::ReversedSequence) = _iscompact(s, s.f)
_iscompact(s::ReversedSequence, f) = iscompact(f)

nzrange(s::ReversedSequence) = _nzrange(s, s.f)

function _nzrange(s::ReversedSequence, f)
    @assert _iscompact(s, f)
    r = nzrange(f)
    -last(r):-first(r)
end

ztransform(s::ReversedSequence, z) = _ztransform(s, z, s.f)
_ztransform(s::ReversedSequence, z, f) = ztransform(f, 1/z)


"""
A DownsampledSequence is the lazy representation of a sequence that is
downsampled by a factor `M`. It is defined by:
`g[k] = f[M*k]`.
"""
immutable DownsampledSequence{T} <: Sequence{T}
    f   ::  Sequence{T}
    M   ::  Int
end

samplefactor(s::DownsampledSequence) = s.M

"Return the sequence downsampled by a factor `M`. By default, `M=2`."
downsample(s::Sequence, M = 2) = DownsampledSequence(s, M)

downsample(s::DownsampledSequence, M = 2) = DownsampledSequence(s.f, s.M*M)

getindex(s::DownsampledSequence, k) = _getindex(s, k, s.f, s.M)
_getindex(s::DownsampledSequence, k, f, M) = f[M*k]

iscompact(s::DownsampledSequence) = _iscompact(s, s.f)
_iscompact(s::DownsampledSequence, f) = iscompact(f)

nzrange(s::DownsampledSequence) = _nzrange(s, s.f, s.M)

function _nzrange(s::DownsampledSequence, f, M)
    @assert _iscompact(s, f)
    r = nzrange(f)
    r1 = first(r)
    r2 = last(r)
    div(r1-1,M)+1 : div(r2-1,M)+1
end

ztransform(s::DownsampledSequence, z) = _ztransform(s, z, s.f, s.M)

function _ztransform(s::DownsampledSequence, z, f, M)
    T = promote_type(eltype(s), eltype(z), eltype(1/z))
    u = zero(T)
    for m = 0:M-1
        u += 1/M * ztransform(f, exp(-im*2*T(pi)*m/M)*z^(1/M))
    end
    u
end


"""
An UpsampledSequence is the lazy representation of a sequence that is
upsampled by a factor `M`, with intermediate zeros. It is defined by:
- `g[k] = f[k/M]`, if `k` is a multiple of `M`, e.g. `k = lM`
- `g[k] = 0`, otherwise.
"""
immutable UpsampledSequence{T} <: Sequence{T}
    f   ::  Sequence{T}
    M   ::  Int
end

samplefactor(s::UpsampledSequence) = s.M

"Return the sequence upsampled by a factor `M`. By default, `M=2`."
upsample(s::Sequence, M = 2) = UpsampledSequence(s, M)

upsample(s::UpsampledSequence, M = 2) = UpsampledSequence(s.f, s.M*M)

getindex(s::UpsampledSequence, k) = _getindex(s, k, s.f, s.M)
_getindex(s::UpsampledSequence, k, f, M) = mod(k, M) == 0 ? f[div(k,M)] : zero(eltype(s))

iscompact(s::UpsampledSequence) = _iscompact(s, s.f)
_iscompact(s::UpsampledSequence, f) = iscompact(f)

nzrange(s::UpsampledSequence) = _nzrange(s, s.f, s.M)

function _nzrange(s::UpsampledSequence, f, M)
    @assert _iscompact(s, f)
    r = nzrange(f)
    r1 = first(r)
    r2 = last(r)
    r1*M : r2*M
end

ztransform(s::UpsampledSequence, z) = _ztransform(s, z, s.f, s.M)
_ztransform(s::UpsampledSequence, z, f, M) = ztransform(f, z^M)


"""
A ModulatedSequence is the lazy representation of a sequence that is modulated:
`g[k] = (-1)^k x[k]`.
"""
immutable ModulatedSequence{T} <: Sequence{T}
    f   ::  Sequence{T}
end

modulate(s::Sequence) = ModulatedSequence(s)
modulate(s::ModulatedSequence) = s.f

getindex(s::ModulatedSequence, k) = _getindex(s, k, s.f)
# We avoid writing (-1)^k, since this throws a DomainError if k is negative
_getindex(s::ModulatedSequence, k, f) = (-1)^isodd(k) * f[k]

iscompact(s::ModulatedSequence) = iscompact(s.f)

nzrange(s::ModulatedSequence) = nzrange(s.f)

ztransform(s::ModulatedSequence, z) = ztransform(s.f, -z)
