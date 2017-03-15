# SignalExtensions.jl

module SignalExtensions

import Base: eltype, getindex, setindex!, eachindex, collect

import Base: &, |, *, transpose, ctranspose, conj

import Base: convert


## We list all functions exported by the package below.

# From sequences.jl:
# - types
export Sequence
# - iterator methods
export firstindex, lastindex, each_nonzero_index
# - generic methods
export moment

# from transforms.jl:
export ztransform, fouriertransform


# From extensions.jl
export PeriodicExtension, ZeroPadding, ConstantPadding, SymmetricExtension

export symmetric_extension_wholepoint_even, symmetric_extension_wholepoint_odd,
    symmetric_extension_halfpoint_even, symmetric_extension_halfpoint_odd


# Definition of abstract Sequence type
include("sequences.jl")

# Fourier and z-transforms
include("transforms.jl")

# The various types of extensions
include("extensions.jl")


end # module
