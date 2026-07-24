# =============================================================================
# Generic Downsampling Engine
# =============================================================================

"""
    downsample(matrix, method::AbstractDownsamplingMethod; kwargs...)

Downsample a binary interaction matrix using the specified algorithm.

Two modes of operation are supported.

- **Single-step downsampling** (`target_co = nothing`)
  Applies a single probabilistic pruning step.

- **Iterative downsampling** (`target_co` specified)
  Repeatedly removes interactions until the requested connectance is
  reached, or no further ecologically valid removals exist.

# Keyword Arguments

- `target_co` : Target connectance.
- `min_spp_prop` : Minimum proportion of active species that must remain.
- `allow_species_loss` : If `false`, no accepted removal may disconnect a
  species.
- `max_iter` : Maximum number of accepted removal attempts.
"""
function downsample(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod;
    target_co::Union{Nothing,Real}=nothing,
    min_spp_prop::Real=0.5,
    allow_species_loss::Bool=true,
    max_iter::Int=500,
)

    if isnothing(target_co)

        return _single_step(matrix, method)

    end

    return _iterative_downsample(
        matrix,
        method;
        target_co,
        min_spp_prop,
        allow_species_loss,
        max_iter,
    )

end


# =============================================================================
# Generic Iterative Engine
# =============================================================================

"""
    _iterative_downsample(...)

Generic connectance-targeted downsampling engine.

This routine contains no method-specific ecological logic. Each
downsampling algorithm simply supplies removal weights via
`_link_removal_weights`.

Ecological constraints (species retention, species loss, etc.) are enforced
centrally by `_choose_removal`.
"""
function _iterative_downsample(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod;
    target_co::Real,
    min_spp_prop::Real,
    allow_species_loss::Bool,
    max_iter::Int,
)

    S = size(matrix, 1)

    current = copy(matrix)

    min_species = ceil(Int, S * min_spp_prop)

    current_co = _connectance(current)

    if current_co <= target_co

        @warn "Initial connectance ($current_co) is already below target ($target_co)."

        return current

    end

    # Best network encountered so far.
    best = copy(current)
    best_co = current_co
    best_diff = abs(current_co - target_co)

    iter = 0

    while current_co > target_co && iter < max_iter

        iter += 1

        # Method-specific removal weights.
        links, weights = _link_removal_weights(current, method)

        isempty(links) && break

        # Select an ecologically valid interaction for removal.
        target_link = _choose_removal(
            current,
            links,
            weights;
            allow_species_loss,
            min_species,
        )

        # No valid removals remain.
        isnothing(target_link) && break

        current[target_link] = false

        current_co = _connectance(current)

        diff = abs(current_co - target_co)

        if diff < best_diff

            best = copy(current)
            best_co = current_co
            best_diff = diff

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

Every downsampling method must implement this function.
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

- `links` contains every realised interaction in the network, and
- `weights` gives the relative probability that each interaction should be
  removed.

Higher weights correspond to a higher probability of removal.

Every iterative downsampling method must implement this function.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    method::AbstractDownsamplingMethod,
)
    throw(MethodError(_link_removal_weights, (matrix, method)))
end


# =============================================================================
# Convenience wrappers
# =============================================================================

# Convert plain Symbol into the Val type for dispatch
downsample(A, method::Symbol; kwargs...) = downsample(A, Val(method); kwargs...)


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

Specify a downsampling algorithm using a symbol.

Examples
--------
```julia
downsample(mat, :powerlaw; y=2.0)

downsample(mat, :niche; sigma_scale=1.0)

downsample(mat, :degree)

downsample(mat, :random)
"""