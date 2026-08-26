using FoodWebTools
using Test
using Random
using Distributions

@testset "Initial body masses" begin
    Random.seed!(66)

    classes = ["small", "medium", "large"]
    bounds = Dict(
        "small" => (1.0, 10.0),
        "medium" => (10.0, 100.0),
        "large" => (100.0, 1000.0)
    )

    dist = LogNormal(2, 1)

    M = FoodWebTools.initial_bodymasses(classes, bounds, dist)

    @test length(M) == 3

    for (m, c) in zip(M, classes)
        lo, hi = bounds[c]
        @test lo <= m <= hi
    end
end

@testset "Log-likelihood" begin
    logM = log.([100.0, 10.0])

    A = [0 1;
         0 0]

    ratio = log(100/10)

    expected = logpdf(Normal(6.1, 5.75), ratio)

    @test FoodWebTools.loglikelihood(logM, A) ≈ expected
end

@testset "Empty food web" begin
    A = zeros(Int, 3, 3)
    logM = rand(3)

    @test FoodWebTools.loglikelihood(logM, A) == 0.0
end

@testset "Proposal bounds" begin
    Random.seed!(66)

    logM = log.([50.0])
    bounds = (10.0, 100.0)

    for _ in 1:1000
        p = FoodWebTools.propose(logM, 1, bounds)
        @test log(bounds[1]) <= p[1] <= log(bounds[2])
    end
end

@testset "Metropolis respects bounds" begin
    Random.seed!(66)

    classes = ["small", "large"]
    bounds = Dict(
        "small" => (1.0, 10.0),
        "large" => (100.0, 1000.0)
    )

    species_bounds = [bounds[c] for c in classes]

    logM = log.([5.0, 200.0])

    A = [0 0;
         1 0]

    for _ in 1:500
        logM = FoodWebTools.metropolis_step(logM, A, species_bounds)
    end

    @test log(bounds["small"][1]) <= logM[1] <= log(bounds["small"][2])
    @test log(bounds["large"][1]) <= logM[2] <= log(bounds["large"][2])
end

@testset "MCMC improves likelihood" begin
    Random.seed!(123)

    classes = fill("medium", 5)
    bounds = Dict("medium" => (1.0, 1000.0))
    dist = LogNormal(2, 1)

    A = [
        0 1 1 0 0;
        0 0 1 0 0;
        0 0 0 1 0;
        0 0 0 0 1;
        0 0 0 0 0
    ]

    M = FoodWebTools.initial_bodymasses(classes, bounds, dist)
    logM = log.(M)
    species_bounds = [bounds[c] for c in classes]

    ll0 = FoodWebTools.loglikelihood(logM, A)

    for _ in 1:5000
        logM = FoodWebTools.metropolis_step(logM, A, species_bounds)
    end

    ll1 = FoodWebTools.loglikelihood(logM, A)

    @test ll1 > ll0
end