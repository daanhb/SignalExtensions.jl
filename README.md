[![Build Status](https://travis-ci.org/daanhb/SignalExtensions.jl.svg?branch=master)](https://travis-ci.org/daanhb/SignalExtensions.jl)
[![Coverage Status](https://coveralls.io/repos/github/daanhb/SignalExtensions.jl/badge.svg)](https://coveralls.io/github/daanhb/SignalExtensions.jl)


# SignalExtensions.jl
=======================

The package provides convenience methods for defining and manipulating bi-infinite digital signals of the form `\{a[k]\}_{k \in \mathbb{Z}}`, where `k` ranges over all integers.

Digital signals can be defined based on finite vectors by periodic extension, symmetric extension (whole and half-point symmetry around either endpoint), zero-padding or constant padding.

Basic support is provided for the Z-transform and Fourier transform of infinite signals.
