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