# test_suite.jl

using SignalExtensions
using Base.Test
using OffsetArrays

include("test_extensions.jl")

function delimit(s::AbstractString)
    println("############")
    println("# ",s)
    println("############")
end

function run_tests()
    delimit("Extensions")
    test_extensions()
end

run_tests()
