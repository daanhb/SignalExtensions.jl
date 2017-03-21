# SignalExtensions.jl

module SignalExtensions

using OffsetArrays

import Base: eltype, getindex, setindex!, eachindex, collect, nzrange

import Base: transpose, ctranspose, conj, reverse

import Base: *


# A collection of useful functions
include("util.jl")
# Definition of abstract Sequence type
include("sequences.jl")
# Definition of z-transform and Fourier transform
include("transforms.jl")
# The various types of extensions
include("extensions.jl")
# Extension types that store a vector
include("extension_sequences.jl")
# Lazy representations of sequence manipulations
include("lazy.jl")
# A number of special sequences
include("special_sequences.jl")

## We list all functions exported by the package below.
## This list is exhaustive and may serve as a reference for the functionality.

# From util.jl
export reverse_indices

# From sequences.jl:
# - types
export Sequence
# - iterator methods
export firstindex, lastindex, each_nonzero_index
# - properties
export iscompact

# From transforms.jl
# - types
export ZTransform, FourierTransform
# - methods
export ztransform, fouriertransform

# From extensions.jl
# - types
export PeriodicExtension, ZeroPadding, ConstantPadding, SymmetricExtension
# - methods
export element, periodic_index

# From extension_sequences.jl
# - types
export ExtensionSequence, PeriodicSequence, CompactSequence, ZeroPaddingSequence,
    ConstantPaddingSequence, SymmetricSequence
# - convenience constructors
export symmetric_extension_wholepoint_even, symmetric_extension_wholepoint_odd,
    symmetric_extension_halfpoint_even, symmetric_extension_halfpoint_odd
# - other methods
export subvector, first_subindex, last_subindex, sublength
export left_parity, right_parity, left_symmetry, right_symmetry

# From lazy.jl
# - types
export ModulatedSequence, DownsampledSequence, UpsampledSequence, ShiftedSequence,
    ReversedSequence, Convolution
# - methods
export downsample, upsample, convolve, shift, modulate, reverse
export sequence_shift, samplefactor

# From special_sequences.jl
export ZeroSequence, DiracSequence

end # module
