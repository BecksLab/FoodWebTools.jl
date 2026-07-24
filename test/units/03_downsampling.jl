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

function active_species(mat)

    rows = vec(sum(mat, dims=2) .> 0)
    cols = vec(sum(mat, dims=1) .> 0)

    sum(rows .|| cols)

end

@testset "Species conservation" begin

    result = downsample(
        TEST_WEB,
        :degree;
        target_co=0.15,
        min_spp_prop=0.8
    )


    @test active_species(result) >= ceil(0.8*5)

end

function connectance(mat)
    sum(mat) / length(mat)
end


@testset "Connectance decreases" begin

    original = connectance(TEST_WEB)

    for method in [
        PowerLaw(2.0),
        Niche(1.0),
        DegreeProduct(),
        RandomSampling()
    ]

        result = downsample(TEST_WEB, method)

        @test connectance(result) <= original

    end

end

@testset "Target connectance" begin

    target = 0.15


    for method in [
        DegreeProduct(),
        RandomSampling(),
        PowerLaw(2.0),
        Niche(1.0)
    ]

        result = downsample(
            TEST_WEB,
            method;
            target_co=target
        )


        co = connectance(result)

        @test co <= connectance(TEST_WEB)

    end

end

# Now we test the ecology

GENERALIST_WEB = Bool[
    1 1 1 1 0
    0 1 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 0
]

@testset "Degree product targets generalists" begin

    removals = 0

    for i in 1:500

        result = downsample(
            GENERALIST_WEB,
            DegreeProduct()
        )

        if result[1, 1] == false
            removals += 1
        end

    end


    @test removals > 0

end

@testset "Power law assigns higher retention to generalists" begin

    WEB = Bool[
        1 1 1 1 1
        0 1 0 0 0
        0 0 1 0 0
        0 0 0 1 0
        0 0 0 0 1
    ]

    probs = FoodWebTools.Downsampling._retention_probabilities(
        WEB,
        PowerLaw(2.0)
    )

    @test probs[1] > probs[2]
    @test probs[1] > probs[3]
    @test probs[1] > probs[4]
    @test probs[1] > probs[5]

end

@testset "Degree product removes highly connected interactions" begin

    DEGREE_WEB = Bool[
        1 1 1 1 1
        1 0 0 0 0
        1 0 0 0 0
        1 0 0 0 0
        1 0 0 0 0
    ]

    Random.seed!(66)

    high_degree_losses = 0
    low_degree_losses = 0


    for i in 1:1000

        result = downsample(
            DEGREE_WEB,
            DegreeProduct()
        )


        if !result[1, 1]
            high_degree_losses += 1
        end

        if !result[2, 1]
            low_degree_losses += 1
        end

    end


    @test high_degree_losses > low_degree_losses

end