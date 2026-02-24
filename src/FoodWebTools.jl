module FoodWebTools

# Dependencies

# export functions

# adbm
include(joinpath("generative_models", "adbm.jl"))
export adbm, adbm_parameters

# l matrix
include(joinpath("generative_models", "lmatrix.jl"))
export lmatrix

# niche
include(joinpath("generative_models", "niche_model.jl"))
export generate_niche_model

# random
include(joinpath("generative_models", "random_model.jl"))
export generate_random_model

# cascade
include(joinpath("generative_models", "cascade_model.jl"))
export generate_cascade_model

# ltm
include(joinpath("generative_models", "ltm.jl"))
export lmatrix

end # module FoodWebTools
