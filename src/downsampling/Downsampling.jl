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
# Method hierarchy
# ============================================================================

"""
Abstract supertype for all network downsampling methods.
"""
abstract type AbstractDownsamplingMethod end

"""
Roopnarine (2006) link-retention model.

`y` controls the scaling of consumer generality.
"""
struct PowerLaw{T<:Real} <: AbstractDownsamplingMethod
    y::T
end

"""
Topological niche-space downsampling.

`sigma_scale` controls the width of consumer niche breadth.
"""
struct Niche{T<:Real} <: AbstractDownsamplingMethod
    sigma_scale::T
end

"""
Degree-product downsampling.

Links involving highly connected consumers and resources are
preferentially removed.
"""
struct DegreeProduct <: AbstractDownsamplingMethod
end

"""
Random link-removal downsampling.

`candidates_per_step` controls how many random candidate removals
are evaluated during iterative downsampling.
"""
struct RandomSampling <: AbstractDownsamplingMethod
    candidates_per_step::Int
end

RandomSampling() = RandomSampling(10)

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