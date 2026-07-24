module Downsampling

using LinearAlgebra
using SpeciesInteractionNetworks
using Statistics
using Random

# ============================================================================
# Public API
# ============================================================================

export downsample

export AbstractDownsamplingMethod
export PowerLaw
export Niche
export DegreeProduct
export RandomSampling

# ============================================================================
# Source files
# ============================================================================

include("methods.jl")
include("utils.jl")
include("engine.jl")

include("power_law.jl")
include("niche.jl")
include("degree_product.jl")
include("random.jl")

end # module