# =============================================================================
# Generic Downsampling Engine
# =============================================================================

"""
    downsample(matrix, method::AbstractDownsamplingMethod; kwargs...)

Downsample a binary interaction matrix using the specified algorithm.

Two modes of operation are supported:

- **Single-step downsampling** (`target_co = nothing`)
    Applies one probabilistic pruning step using the supplied method.

- **Iterative downsampling** (`target_co` specified)
    Repeatedly removes links until the requested connectance is reached,
    or no further valid removals are possible.
"""
function downsample(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod;
    target_co::Union{Nothing,Real}=nothing,
    min_spp_prop::Real=0.5,
    max_iter::Int=500,
)

    # Single-step downsampling
    if isnothing(target_co)
        return _single_step(matrix, method)
    end

    # Iterative connectance-targeted downsampling
    return _iterative_downsample(
        matrix,
        method;
        target_co,
        min_spp_prop,
        max_iter,
    )

end


# =============================================================================
# Generic Iterative Engine
# =============================================================================

"""
    _iterative_downsample(...)

Generic iterative downsampling routine.

This function is completely agnostic to the downsampling algorithm being
used. Algorithm-specific behaviour is supplied through multiple dispatch
via `_select_candidate()`.

Each iteration proceeds as follows:

1. Ask the algorithm which link should be removed.
2. Test whether removing that link violates species-retention constraints.
3. Accept or reject the removal.
4. Track the network closest to the requested connectance.
"""
function _iterative_downsample(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod;
    target_co::Real,
    min_spp_prop::Real,
    max_iter::Int,
)

    S = size(matrix, 1)

    current = copy(matrix)

    min_species = ceil(Int, S * min_spp_prop)

    active_species, current_co = _get_downsample_metrics(current)

    # Nothing to do if the network is already sparse enough.
    if current_co <= target_co
        @warn "Initial connectance ($current_co) is already below target ($target_co)."
        return current
    end

    # Track the closest network encountered.
    best = copy(current)
    best_co = current_co
    best_diff = abs(current_co - target_co)

    iter = 0

    while current_co > target_co && iter < max_iter

        iter += 1

        # Ask the algorithm which interaction to attempt removing.
        #
        # Returning `nothing` signals that no valid candidates remain.
        target_link = _select_candidate(current, method)

        isnothing(target_link) && break

        temp = copy(current)
        temp[target_link] = false

        active_species, next_co = _get_downsample_metrics(temp)

        # Reject removals that eliminate too many species.
        if active_species < min_species
            continue
        end

        # Accept the removal.
        current = temp
        current_co = next_co

        # Update the best network encountered.
        diff = abs(current_co - target_co)

        if diff < best_diff
            best = copy(current)
            best_diff = diff
            best_co = current_co
        end

    end

    if iter == max_iter && current_co > target_co
        @warn "Reached maximum iterations before target connectance was achieved. Closest connectance = $best_co."
    end

    return best

end


# =============================================================================
# Extension Interface
# =============================================================================

"""
    _single_step(matrix, method)

Perform a single probabilistic downsampling step.

Every downsampling algorithm must implement this method.
"""
function _single_step(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod,
)
    throw(MethodError(_single_step, (matrix, method)))
end


"""
    _link_removal_weights(matrix, method)

Return

    links, weights

where

- `links` is the collection of existing interactions.
- `weights` gives each interaction's relative probability of removal.

Most algorithms only need to implement this function.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod,
)
    throw(MethodError(_link_removal_weights, (matrix, method)))
end


"""
    _select_candidate(matrix, method)

Select the interaction that should be tested for removal.

The default implementation performs weighted random sampling using
the weights returned by `_link_removal_weights()`.

Algorithms with more sophisticated candidate-selection strategies
may overload this function.
"""
function _select_candidate(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod,
)

    links, weights = _link_removal_weights(matrix, method)

    isempty(links) && return nothing

    probabilities = weights ./ sum(weights)

    return links[_rand_categorical(probabilities)]

end


# =============================================================================
# Convenience Wrappers
# =============================================================================

downsample(mat, ::Val{:powerlaw}; y, kwargs...) =
    downsample(mat, PowerLaw(y); kwargs...)

downsample(mat, ::Val{:niche}; sigma_scale, kwargs...) =
    downsample(mat, Niche(sigma_scale); kwargs...)

downsample(mat, ::Val{:degree}; kwargs...) =
    downsample(mat, DegreeProduct(); kwargs...)

downsample(mat, ::Val{:random}; kwargs...) =
    downsample(mat, RandomSampling(); kwargs...)

"""
    downsample(matrix, method::Symbol; kwargs...)

Convenience wrapper allowing algorithms to be specified by name.

Examples
--------
```julia
downsample(mat, :powerlaw; y=2.0)

downsample(mat, :niche; sigma_scale=0.8)

downsample(mat, :random)
```
"""
function downsample(mat, method::Symbol; kwargs...)

    return downsample(mat, Val(method); kwargs...)

end