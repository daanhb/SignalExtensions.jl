# test_lazy.jl

function test_lazy()
    @testset "$(rpad("Shifted sequences",80))" begin
        test_shift() end
    @testset "$(rpad("Reversed sequences",80))" begin
        test_reverse() end
    @testset "$(rpad("Downsampled sequences",80))" begin
        test_downsample() end
    @testset "$(rpad("Upsampled sequences",80))" begin
        test_upsample() end
    @testset "$(rpad("Modulated sequences",80))" begin
        test_modulate() end
    @testset "$(rpad("Convolutions",80))" begin
        test_convolution() end
end

function test_shift()
    a = CompactSequence([1,2,3,4,5,6])
    test_shifted_sequence(shift(a, 2), a)
    test_shifted_sequence(ShiftedSequence(a, 2), a)
end

function test_shifted_sequence(s, f)
    n = sequence_shift(s)

    @test s[4] == f[4-n]

    @test iscompact(s) == iscompact(f)

    if iscompact(f)
        @test first(nzrange(s)) == first(nzrange(f))+n
        @test last(nzrange(s)) == last(nzrange(f))+n
    end

    # Check that an additional shift is added to the existing shift
    if typeof(s) <: ShiftedSequence
        @test sequence_shift(shift(s, 4)) == n+4
    end
end

function test_reverse()
    a = CompactSequence([1,2,3,4,5,6])
    test_reversed_sequence(reverse(a), a)
    test_reversed_sequence(ReversedSequence(a), a)
end

function test_reversed_sequence(s, f)
    @test s[4] == f[-4]
    @test iscompact(s) == iscompact(f)

    if iscompact(f)
        @test first(nzrange(s)) == -last(nzrange(f))
        @test last(nzrange(s)) == -first(nzrange(f))
    end

    if typeof(s) <: ReversedSequence
        @test reverse(s) == f
    end
end

function test_downsample()
    a = CompactSequence([1,2,3,4,5,6])
    test_downsampled_sequence(downsample(a, 2), a, 2)
    test_downsampled_sequence(DownsampledSequence(a, 2), a, 2)
end

function test_downsampled_sequence(s, f, M)
    @test samplefactor(s) == M
    for k in (-2, 0, 1, 4)
        @test s[k] == f[M*k]
    end
end

function test_upsample()
    a = CompactSequence([1,2,3,4,5,6])
    test_upsampled_sequence(upsample(a, 2), a, 2)
    test_upsampled_sequence(UpsampledSequence(a, 2), a, 2)
end

function test_upsampled_sequence(s, f, M)
    @test samplefactor(s) == M
    for k in (-2, 0, 1, 4)
        @test s[M*k] == f[k]
    end
    for m in 1:M-1
        @test s[m] == zero(eltype(s))
    end
end

function test_modulate()
    a = CompactSequence([1.0,2.0,3.0,4.0,5.0,6.0])
    test_modulated_sequence(modulate(a), a)
    test_modulated_sequence(ModulatedSequence(a), a)
end

function test_modulated_sequence(s, f)
    for k in -2:2
        @test s[k] == (-one(eltype(s)) )^k * f[k]
    end

    @test iscompact(s) == iscompact(f)

    if iscompact(f)
        @test first(nzrange(s)) == first(nzrange(f))
        @test last(nzrange(s)) == last(nzrange(f))
    end
end

function test_convolution()
    a = CompactSequence([1.0,2.0,3.0,4.0,5.0,6.0])
    b = CompactSequence(OffsetArray([4.1,2.2,3.3], -2:0))
    @test a*a == convolve(a, a)
    test_convolved_sequences(a*a, a, a)
    test_convolved_sequences(a*b, a, b)
end

function test_convolved_sequences(s, f, g)
    n = 4
    z = zero(eltype(f))
    @test s[n] == sum([f[k]*g[n-k] for k in nzrange(f)])

    @test iscompact(s) == iscompact(f) && iscompact(g)
end
