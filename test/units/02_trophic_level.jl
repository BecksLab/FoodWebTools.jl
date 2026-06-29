using Test
using FoodWebTools

@testset "trophic_level" begin

    @testset "Producer → Consumer" begin
        # Species 2 eats species 1
        A = Bool[
            0 0;
            1 0
        ]

        tl = trophic_level(A)

        @test tl ≈ [1.0, 2.0]
    end

    @testset "Three-species chain" begin
        # 3 → 2 → 1
        A = Bool[
            0 0 0;
            1 0 0;
            0 1 0
        ]

        tl = trophic_level(A)

        @test tl ≈ [1.0, 2.0, 3.0]
    end

    @testset "Omnivore" begin
        # Species 3 eats species 1 and 2
        #
        # 1 basal
        # 2 eats 1
        # 3 eats 1 and 2
        #
        # TL(3) = 1 + (1 + 2)/2 = 2.5

        A = Bool[
            0 0 0;
            1 0 0;
            1 1 0
        ]

        tl = trophic_level(A)

        @test tl ≈ [1.0, 2.0, 2.5]
    end

    @testset "Named species" begin
        A = Bool[
            0 0;
            1 0
        ]

        names = ["Plant", "Herbivore"]

        tl = trophic_level(A; species=names)

        @test tl isa Dict
        @test tl["Plant"] ≈ 1.0
        @test tl["Herbivore"] ≈ 2.0
    end

    @testset "Output length" begin
        A = rand(Bool, 10, 10)

        tl = trophic_level(A)

        @test length(tl) == 10
    end

end