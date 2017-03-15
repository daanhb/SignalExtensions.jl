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

function test_vectors()
    test_vectors = Array(Any, 0)
    # A regular array
    push!(test_vectors, [1,2,3])
    # A LinSpace object
    push!(test_vectors, [1,2,3])
    push!(test_vectors, OffsetArray([1,2,3,4,5], -1:3))
    test_vectors
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
    for i in eachindex(a)
        @test p[i] == a[i]
    end
    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == a[lastindex]
    @test p[firstindex-2] == a[lastindex-1]
    @test p[lastindex+1] == a[firstindex]
    @test p[lastindex+2] == a[firstindex+1]

    some_index = 100
    periodic_index = mod(some_index-firstindex, n) + firstindex
    @test p[some_index] == a[periodic_index]
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
    for i in eachindex(a)
        @test p[i] == a[i]
    end
    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == 0
    @test p[firstindex-2] == 0
    @test p[lastindex+1] == 0
    @test p[lastindex+2] == 0

    # Check that the zerp has the same type as a
    @test eltype(p[firstindex-1]) == eltype(a)
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
    for i in eachindex(a)
        @test p[i] == a[i]
    end
    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1
    @test p[firstindex-1] == c
    @test p[firstindex-2] == c
    @test p[lastindex+1] == c
    @test p[lastindex+2] == c

    # Check that the constant has the same type as a
    @test eltype(p[firstindex-1]) == eltype(a)
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
    for i in eachindex(a)
        @test p[i] == a[i]
    end
    firstindex = first(linearindices(a))
    lastindex = last(linearindices(a))
    n = lastindex-firstindex+1

    p1 = symmetric_extension_wholepoint_even(a)
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
end
