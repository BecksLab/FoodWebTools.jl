using FoodWebTools
using Test

const TEST_WEB = Bool[
    0 1 1 0 0
    0 0 1 1 0
    0 0 0 1 1
    0 0 0 0 1
    1 0 0 0 0
]


@testset "Downsampling API" begin

    matrix = TEST_WEB


    @testset "Single-step methods return matrices" begin

        results = [
            downsample(matrix, :random),
            downsample(matrix, :degree),
            downsample(matrix, :powerlaw; y=2.0),
            downsample(matrix, :niche; sigma_scale=1.0),
        ]

        for i in eachindex(results)

            @test results[i] isa AbstractMatrix{Bool}
            @test size(results[i]) == size(matrix)

        end

    end


    @testset "No new links are created" begin

        for result in [
            downsample(matrix, :random),
            downsample(matrix, :degree),
            downsample(matrix, :powerlaw; y=2.0),
            downsample(matrix, :niche; sigma_scale=1.0),
        ]

            @test all(result .<= matrix)

        end

    end


end

@testset "Target connectance downsampling" begin

    target = 0.10

    for method in [
        :random,
        :degree,
        :powerlaw,
        :niche,
    ]

        result =
            if method == :powerlaw
                downsample(TEST_WEB, method;
                    y=2.0,
                    target_co=target
                )

            elseif method == :niche
                downsample(TEST_WEB, method;
                    sigma_scale=1.0,
                    target_co=target
                )

            else
                downsample(TEST_WEB, method;
                    target_co=target
                )
            end


        @test size(result) == size(TEST_WEB)

        # Should never exceed original links
        @test all(result .<= TEST_WEB)


        # Connectance should have decreased
        original_links = sum(TEST_WEB)
        new_links = sum(result)

        @test new_links <= original_links

    end

end

@testset "Species retention constraint" begin

    result = downsample(
        TEST_WEB,
        :degree;
        target_co=0.05,
        min_spp_prop=0.8
    )


    active_species =
        sum(
            vec(sum(result, dims=1) .> 0) .||
            vec(sum(result, dims=2) .> 0)
        )


    @test active_species >= ceil(Int, 5*0.8)

end

using Random

@testset "Random reproducibility" begin

    Random.seed!(66)

    a = downsample(TEST_WEB, :random)

    Random.seed!(66)

    b = downsample(TEST_WEB, :random)

    @test a == b

end