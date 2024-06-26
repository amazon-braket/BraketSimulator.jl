using Test, BraketSimulator

@testset "utils" begin
    for ix = 0:(2^12-1)
        @test BraketSimulator._index_to_endian_bits(ix, 12) ==
              reverse(digits(ix, base = 2, pad = 12))
    end
end
