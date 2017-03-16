# util.jl

"""
Return a vector with reversed indices, if `p = reverse_indices(a)` then
`p[i] = a[-i]`. The result is an OffsetArray.
"""
function reverse_indices(a::AbstractVector)
    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    OffsetArray([a[i] for i in lastindex:-1:firstindex], -lastindex:-firstindex)
end
