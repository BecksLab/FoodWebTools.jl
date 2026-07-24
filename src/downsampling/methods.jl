# =============================================================================
# Downsampling Methods
# =============================================================================

"""
Abstract supertype for all network downsampling methods.

New downsampling algorithms should subtype this and implement

- `_single_step`
- `_link_removal_weights`

to integrate with the generic downsampling engine.
"""
abstract type AbstractDownsamplingMethod end

# =============================================================================
# PowerLaw - Roopnarine (2006)
# =============================================================================

"""
    PowerLaw(y)

Downsampling based on the scaling model of Roopnarine (2006).

The parameter `y` controls how strongly consumer generality influences
link retention.
"""
struct PowerLaw{T<:Real} <: AbstractDownsamplingMethod
    y::T
end

# =============================================================================
# Topological niche
# =============================================================================

"""
    Niche(sigma_scale)

Downsampling based on topological niche positions derived from the
leading singular vectors of the interaction matrix.

`sigma_scale` controls the effective niche breadth.
"""
struct Niche{T<:Real} <: AbstractDownsamplingMethod
    sigma_scale::T
end

# =============================================================================
# Degree-product
# =============================================================================

"""
    DegreeProduct()

Downsampling where links connecting highly connected consumers and
resources are preferentially removed.
"""
struct DegreeProduct <: AbstractDownsamplingMethod
end

# =============================================================================
# Random
# =============================================================================

"""
    RandomSampling([candidates_per_step])

Random link removal.

For iterative downsampling, `candidates_per_step` random links are
evaluated at each iteration and the best candidate is selected.
"""
struct RandomSampling <: AbstractDownsamplingMethod
    candidates_per_step::Int
end