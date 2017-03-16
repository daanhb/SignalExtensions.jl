# test_extensions.jl

function test_extensions()
    @testset "$(rpad("Periodic extension",80))" begin
        test_periodic_extensions() end
    @testset "$(rpad("Zero padding",80))" begin
        test_zeropadding_extensions() end
    @testset "$(rpad("Constant padding",80))" begin
        test_constantpadding_extensions() end
    @testset "$(rpad("Symmetric extension",80))" begin
        test_symmetric_extensions() end
end

# A collection of vectors of different types for use in the tests
function test_vectors()
    test_vectors = Array(Any, 0)
    # A regular array
    push!(test_vectors, [1,2,3])
    push!(test_vectors, [1.0+1.0im, 2.0+3.0im, 4.2-2im])
    # A LinSpace object
    push!(test_vectors, [1,2,3])
    push!(test_vectors, OffsetArray([1,2,3,4,5], -1:3))
    test_vectors
end

# Test assignment: we only test a true Vector, since one can not assign to
test_assignment(s::ExtensionSequence) = _test_assignment(s, subvector(s))

_test_assignment(s::ExtensionSequence, a::AbstractVector) = none
_test_assignment(s::ExtensionSequence, a::OffsetArray) = do_test_assignment(s)
_test_assignment(s::ExtensionSequence, a::Vector) = do_test_assignment(s)

function do_test_assignment(s::ExtensionSequence)
    val = one(eltype(s))
    idx = first_subindex(s)+1
    s[idx] = val
    @test s[idx] == val
end

function test_periodic_extensions()
    vectors = test_vectors()
    for vector in vectors
        test_periodic_extension(vector)
    end
end

function test_periodic_extension(a)
    p = PeriodicExtension(a)

    @test eltype(p) == eltype(a)

    @test !iscompact(p)

    # Equality of the embedded elements
    for i in eachindex(a)
        @test p[i] == a[i]
    end

    test_assignment(p)

    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == a[lastindex]
    @test p[firstindex-2] == a[lastindex-1]
    @test p[lastindex+1] == a[firstindex]
    @test p[lastindex+2] == a[firstindex+1]

    large_index = 100
    periodic_index = mod(large_index-firstindex, n) + firstindex
    @test p[large_index] == a[periodic_index]

    # Check time-reversal and (c)transpose
    test_index = firstindex + 1
    @test p'[-test_index] == conj(p[test_index])
    @test p.'[-test_index] == p[test_index]
end

function test_zeropadding_extensions()
    vectors = test_vectors()
    for vector in vectors
        test_zeropadding_extension(vector)
    end
end

function test_zeropadding_extension(a)
    p = ZeroPadding(a)

    @test eltype(p) == eltype(a)

    @test iscompact(p)

    # Equality of the embedded elements
    for i in eachindex(a)
        @test p[i] == a[i]
    end

    test_assignment(p)

    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == 0
    @test p[firstindex-2] == 0
    @test p[lastindex+1] == 0
    @test p[lastindex+2] == 0

    # Check that the zero has the same type as a
    @test eltype(p[firstindex-1]) == eltype(a)

    # Check time-reversal and (c)transpose
    test_index = firstindex + 1
    @test p'[-test_index] == conj(p[test_index])
    @test p.'[-test_index] == p[test_index]
end

function test_constantpadding_extensions()
    vectors = test_vectors()
    for vector in vectors
        test_constantpadding_extension(vector, 1)
    end
end

function test_constantpadding_extension(a, c)
    p = ConstantPadding(a, c)

    @test eltype(p) == eltype(a)

    @test !iscompact(p)

    # Equality of the embedded elements
    for i in eachindex(a)
        @test p[i] == a[i]
    end

    test_assignment(p)

    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == c
    @test p[firstindex-2] == c
    @test p[lastindex+1] == c
    @test p[lastindex+2] == c

    # Check that the constant has the same type as a
    @test eltype(p[firstindex-1]) == eltype(a)

    # Check time-reversal and (c)transpose
    test_index = firstindex + 1
    @test p'[-test_index] == conj(p[test_index])
    @test p.'[-test_index] == p[test_index]
end

function test_symmetric_extensions()
    vectors = test_vectors()
    for vector in vectors
        test_symmetric_extension(vector)
    end
end

function test_symmetric_extension(a)
    p = SymmetricExtension(a)

    @test eltype(p) == eltype(a)

    @test !iscompact(p)

    # Equality of the embedded elements
    for i in eachindex(a)
        @test p[i] == a[i]
    end

    test_assignment(p)

    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1

    p1 = symmetric_extension_wholepoint_even(a)
    @test left_parity(p1) == :even
    @test right_parity(p1) == :even
    @test left_symmetry(p1) == :wp
    @test right_symmetry(p1) == :wp
    @test p1[lastindex+1] == a[lastindex]
    @test p1[lastindex+2] == a[lastindex-1]
    @test p1[firstindex-1] == a[firstindex]
    @test p1[firstindex-2] == a[firstindex+1]
    p2 = symmetric_extension_wholepoint_odd(a)
    @test p2[lastindex+1] == -a[lastindex]
    @test p2[lastindex+2] == -a[lastindex-1]
    @test p2[firstindex-1] == -a[firstindex]
    @test p2[firstindex-2] == -a[firstindex+1]
    p3 = symmetric_extension_halfpoint_even(a)
    @test p3[lastindex+1] == a[lastindex-1]
    @test p3[firstindex-1] == a[firstindex+1]
    p4 = symmetric_extension_halfpoint_odd(a)
    @test p4[lastindex+1] == -a[lastindex-1]
    @test p4[firstindex-1] == -a[firstindex+1]

    # A few mixed cases
    p5 = SymmetricExtension{:wp,:hp,:even,:odd,typeof(a),eltype(a)}(a)
    @test p5[lastindex+1] == -a[lastindex-1]
    @test p5[firstindex-1] == a[firstindex]

    p6 = SymmetricExtension{:hp,:wp,:odd,:even,typeof(a),eltype(a)}(a)
    @test p6[lastindex+1] ==  a[lastindex]
    @test p6[firstindex-1] == -a[firstindex+1]

    # Check time-reversal and (c)transpose
    test_index = firstindex + 1
    @test p'[-test_index] == conj(p[test_index])
    @test p.'[-test_index] == p[test_index]
end
