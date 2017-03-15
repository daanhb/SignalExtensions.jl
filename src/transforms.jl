# transforms.jl

"""
Compute the Z transform of a sequence, defined as `S(z) = \sum_k s_k z^{-k}`.
"""
function ztransform(s::Sequence, z)
    T = promote_type(eltype(s), eltype(z), eltype(1/z))
    S = zero(T)
    for k in each_nonzero_index(s)
        S += s[k] * z^(-k)
    end
    S
end

"""
The Fourier transform of a sequence is defined as `S(ω) = \sum_k s_k e^{-i ω k}`.
It corresponds to the Z transform with `z = e^{i ω}`. The Fourier transform is a
`2π`-periodic function of `ω`.
"""
fouriertransform(s::Sequence, ω) = ztransform(s, exp(im*ω))


"""
The ZTransform type is a lazy wrapper for the `ztransform` function.
"""
immutable ZTransform{S <: Sequence}
    sequence     ::  S
end

sequence(zt::ZTransform) = zt.sequence

(zt::ZTransform)(z) = ztransform(sequence(zt), z)

"""
The FourierTransform type is a lazy wrapper for the `fouriertransform` function.
"""
immutable FourierTransform{S <: Sequence}
    seq     ::  S
end

sequence(ft::FourierTransform) = ft.sequence

(ft::FourierTransform)(ω) = fouriertransform(sequence(ft), ω)
