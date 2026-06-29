using Test
using FoodWebTools

@testset "Extinctions" begin

     @testset "All Good" begin
        include("units/00_allgood.jl")
    end

    @testset "Intervality" begin
        include("units/01_intervality.jl")
    end

    @testset "Trophic level" begin
        include("units/02_trophic_level.jl")
    end

end