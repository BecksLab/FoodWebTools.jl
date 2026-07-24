using FoodWebTools
using Test
using Random
using Statistics

# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

const TEST_WEB = Bool[
    0 1 1 0 0
    0 0 1 1 0
    0 0 0 1 1
    0 0 0 0 1
    1 0 0 0 0
]

active_species(mat) =
    sum(
        vec(sum(mat, dims=2) .> 0) .||
        vec(sum(mat, dims=1) .> 0)
    )

connectance(mat) = sum(mat) / length(mat)

# -------------------------------------------------------------------------
# Generic API
# -------------------------------------------------------------------------

@testset "Downsampling API" begin

    for method in (
        :random,
        :degree,
        :powerlaw,
        :niche,
    )

        result =
            method == :powerlaw ? downsample(TEST_WEB, method; y=2.0) :
            method == :niche    ? downsample(TEST_WEB, method; sigma_scale=1.0) :
                                  downsample(TEST_WEB, method)

        @test size(result) == size(TEST_WEB)
        @test result isa AbstractMatrix{Bool}
        @test all(result .<= TEST_WEB)

    end

end

# -------------------------------------------------------------------------
# Target connectance
# -------------------------------------------------------------------------

@testset "Target connectance" begin

    target = 0.15

    for method in (
        :random,
        :degree,
        :powerlaw,
        :niche,
    )

        result =
            method == :powerlaw ? downsample(TEST_WEB, method; y=2.0, target_co=target) :
            method == :niche    ? downsample(TEST_WEB, method; sigma_scale=1.0, target_co=target) :
                                  downsample(TEST_WEB, method; target_co=target)

        @test connectance(result) <= connectance(TEST_WEB)
        @test all(result .<= TEST_WEB)

    end

end

# -------------------------------------------------------------------------
# Species retention
# -------------------------------------------------------------------------

@testset "Minimum species proportion respected" begin

    result = downsample(
        TEST_WEB,
        :degree;
        target_co=0.05,
        min_spp_prop=0.8,
    )

    @test active_species(result) >= 4

end

# -------------------------------------------------------------------------
# New guard rail
# -------------------------------------------------------------------------

@testset "Species loss protection" begin

    WEB = Bool[
        1 0 0
        0 1 0
        0 0 1
    ]

    protected = downsample(
        WEB,
        :degree;
        target_co=0.0,
        allow_species_loss=false,
    )

    unprotected = downsample(
        WEB,
        :degree;
        target_co=0.0,
        allow_species_loss=true,
    )

    @test protected == WEB
    @test sum(unprotected) <= sum(WEB)

end

# -------------------------------------------------------------------------
# Random reproducibility
# -------------------------------------------------------------------------

@testset "Random seed reproducibility" begin

    for method in (
        :random,
        :degree,
        :powerlaw,
        :niche,
    )

        Random.seed!(66)

        a =
            method == :powerlaw ? downsample(TEST_WEB, method; y=2.0) :
            method == :niche    ? downsample(TEST_WEB, method; sigma_scale=1.0) :
                                  downsample(TEST_WEB, method)

        Random.seed!(66)

        b =
            method == :powerlaw ? downsample(TEST_WEB, method; y=2.0) :
            method == :niche    ? downsample(TEST_WEB, method; sigma_scale=1.0) :
                                  downsample(TEST_WEB, method)

        @test a == b

    end

end

# -------------------------------------------------------------------------
# Empty network
# -------------------------------------------------------------------------

@testset "Empty networks remain empty" begin

    WEB = falses(5,5)

    for method in (
        :random,
        :degree,
        :powerlaw,
        :niche,
    )

        result =
            method == :powerlaw ? downsample(WEB, method; y=2.0) :
            method == :niche    ? downsample(WEB, method; sigma_scale=1.0) :
                                  downsample(WEB, method)

        @test result == WEB

    end

end

# -------------------------------------------------------------------------
# Complete network
# -------------------------------------------------------------------------

@testset "Complete networks only lose links" begin

    WEB = trues(5,5)

    for method in (
        :random,
        :degree,
        :powerlaw,
        :niche,
    )

        result =
            method == :powerlaw ? downsample(WEB, method; y=2.0) :
            method == :niche    ? downsample(WEB, method; sigma_scale=1.0) :
                                  downsample(WEB, method)

        @test all(result .<= WEB)

    end

end

# -------------------------------------------------------------------------
# Ecology: degree product
# -------------------------------------------------------------------------

@testset "Degree product preferentially removes hub interactions" begin

    WEB = Bool[
        1 1 1 1 1
        1 0 0 0 0
        1 0 0 0 0
        1 0 0 0 0
        1 0 0 0 0
    ]

    hub_losses = 0
    peripheral_losses = 0

    for i in 1:1000

        result = downsample(WEB, :degree)

        hub_losses += !result[1,1]
        peripheral_losses += !result[2,1]

    end

    @test hub_losses > peripheral_losses

end

# -------------------------------------------------------------------------
# Ecology: power law
# -------------------------------------------------------------------------

@testset "Power law retains generalists" begin

    WEB = Bool[
        1 1 1 1 1
        0 1 0 0 0
        0 0 1 0 0
        0 0 0 1 0
        0 0 0 0 1
    ]

    generalist = Float64[]
    specialist = Float64[]

    for i in 1:500

        result = downsample(WEB, :powerlaw; y=2.0)

        push!(generalist, sum(result[1,:]))

        push!(
            specialist,
            mean(sum.(eachrow(result)[2:5]))
        )

    end

    @test mean(generalist) > mean(specialist)

end

# -------------------------------------------------------------------------
# Already below target
# -------------------------------------------------------------------------

@testset "Already below target connectance" begin

    result = downsample(
        TEST_WEB,
        :degree;
        target_co=1.0,
    )

    @test result == TEST_WEB

end