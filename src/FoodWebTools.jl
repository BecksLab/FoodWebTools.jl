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

# Downsampling
include(joinpath("downsampling", "Downsampling.jl")) # note calls exports already

# Public API
export adbm, adbm_parameters
export lmatrix
export niche_model
export random_model
export cascade_model
export ltm

# Network metrics
export intervality
export trophic_level

using .Downsampling

export downsample
export AbstractDownsamplingMethod
export PowerLaw
export Niche
export DegreeProduct
export RandomSampling

end