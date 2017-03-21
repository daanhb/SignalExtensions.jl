# test_suite.jl

using SignalExtensions
using Base.Test
using OffsetArrays

include("test_util.jl")
include("test_sequences.jl")
include("test_extensions.jl")
include("test_lazy.jl")

function delimit(s::AbstractString)
    println("############")
    println("# ",s)
    println("############")
end

function run_tests()
    delimit("Utility functions")
    test_util()

    delimit("Sequences")
    test_sequences()

    delimit("Extensions")
    test_extensions()

    delimit("Lazy sequences")
    test_lazy()
end

run_tests()
