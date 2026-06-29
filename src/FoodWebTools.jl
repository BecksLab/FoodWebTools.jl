module FoodWebTools

using LinearAlgebra

# Generative models
include(joinpath("generators", "adbm.jl"))
include(joinpath("generators", "lmatrix.jl"))
include(joinpath("generators", "niche.jl"))
include(joinpath("generators", "random.jl"))
include(joinpath("generators", "cascade.jl"))
include(joinpath("generators", "ltm.jl"))

# Network metrics
include(joinpath("metrics", "intervality.jl"))
include(joinpath("metrics", "trophic_level.jl"))

# Public API
export adbm, adbm_parameters
export lmatrix
export niche_model
export random_model
export cascade_model
export ltm

export intervality
export trophic_level

end