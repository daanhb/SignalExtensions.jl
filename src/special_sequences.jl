# special_sequences.jl
# A collection of special sequences.

"The ZeroSequence is a sequence of all zeros."
immutable ZeroSequence{T} <: Sequence{T}
end

ZeroSequence() = ZeroSequence{Float64}()

getindex(s::ZeroSequence, k) = zero(eltype(s))

ztransform(s::ZeroSequence, z) = zero(eltype(s))


"The DiracSequence is a sequence of all zeros, except for `h[0] = 1`."
immutable DiracSequence{T} <: Sequence{T}
end

DiracSequence() = DiracSequence{Float64}()

getindex(s::DiracSequence, k) = k == 0 ? one(eltype(s)) : zero(eltype(s))

ztransform(s::DiracSequence, z) = one(eltype(s))
