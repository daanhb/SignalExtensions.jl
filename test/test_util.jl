# test_util.jl

function test_util()
    @testset "$(rpad("Index reversal",80))" begin
        test_reverse([1,2,3])
        test_reverse(OffsetArray([1,2,3], -2:0))
    end
end

function test_reverse(a)
    b = reverse_indices(a)
    for i in eachindex(a)
        @test a[i] == b[-i]
    end
end
