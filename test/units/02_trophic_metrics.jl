@testset "trophic_level" begin

    @testset "basal species" begin
        A = falses(3,3)

        tl = trophic_level(A)

        @test all(tl .≈ 1)
    end


    @testset "linear food chain" begin
        A = Bool[
            false true false
            false false true
            false false false
        ]

        tl = trophic_level(A)

        @test tl[1] ≈ 3
        @test tl[2] ≈ 2
        @test tl[3] ≈ 1
    end

end


@testset "trophic_coherence" begin

    @testset "perfect chain" begin
        A = Bool[
            false true false
            false false true
            false false false
        ]

        @test trophic_coherence(A) ≈ 0
    end


    @testset "shortcut predation" begin
        A = Bool[
            false true true
            false false true
            false false false
        ]

        @test trophic_coherence(A) > 0
    end


    @testset "empty network" begin
        @test isnan(trophic_coherence(falses(4,4)))
    end

end