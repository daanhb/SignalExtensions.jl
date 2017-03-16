# SignalExtensions.jl

module SignalExtensions

using OffsetArrays

import Base: eltype, getindex, setindex!, eachindex, collect

import Base: &, |, *, transpose, ctranspose, conj, reverse

import Base: convert


# A collection of useful functions
include("util.jl")
# Definition of abstract Sequence type
include("sequences.jl")
# The various types of extensions
include("extensions.jl")


## We list all functions exported by the package below.
## This list is exhaustive and may serve as a reference.

# From util.jl
export reverse_indices

# From sequences.jl:
# - types
export Sequence
# - iterator methods
export firstindex, lastindex, each_nonzero_index
# - properties
export iscompact


# From extensions.jl
# - types
export PeriodicExtension, ZeroPadding, ConstantPadding, SymmetricExtension
# - convenience constructors
export symmetric_extension_wholepoint_even, symmetric_extension_wholepoint_odd,
    symmetric_extension_halfpoint_even, symmetric_extension_halfpoint_odd


end # module
