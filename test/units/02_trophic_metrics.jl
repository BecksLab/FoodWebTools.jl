@testset "Trophic Level Tests" begin

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


    @testset "cannibalism excluded by default" begin
        A = Bool[
            true  true false
            false false true
            false false false
        ]

        tl = trophic_level(A)

        # Species 1 consumes species 2 and itself.
        # The self-link should be ignored, leaving a two-step chain.
        @test tl[1] ≈ 3
        @test tl[2] ≈ 2
        @test tl[3] ≈ 1
    end


    @testset "cannibalism can be retained" begin
        A = Bool[
            true  true false
            false false true
            false false false
        ]

        tl = trophic_level(A; exclude_cannibalism=false)

        # Retaining self-links changes the trophic calculation.
        @test tl[1] != 3
    end


    @testset "cannibalistic-only species treated as basal" begin
        A = Bool[
            true false false
            false false true
            false false false
        ]

        tl = trophic_level(A)

        # Species 1 only consumes itself. Removing the self-link leaves it
        # without prey, so it is assigned trophic level 1.
        @test tl[1] ≈ 1
    end


    @testset "input adjacency matrix is unchanged" begin
        A = Bool[
            true  true false
            false false true
            false false false
        ]

        A_original = copy(A)

        trophic_level(A)

        @test A == A_original
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