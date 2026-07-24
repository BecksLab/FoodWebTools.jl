using Test
using FoodWebTools

@testset "FoodWeb tools" begin

     @testset "All Good" begin
        include("units/00_allgood.jl")
    end

    @testset "Trophic metrics" begin
        include("units/02_trophic_metrics.jl")
    end

    @testset "Downsampling" begin
        include("units/03_downsampling.jl")
    end

    @testset "Network metrics" begin
        include("units/04_network_metrics.jl")
    end

end