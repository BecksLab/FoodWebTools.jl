# =============================================================================
# Shared Utility Functions
# =============================================================================

"""
    _get_downsample_metrics(mat)

Calculate active species count and network connectance.

Returns:

    (active_species, connectance)

where connectance is:

    L / S²

with `L` the number of realised interactions.
"""
function _get_downsample_metrics(
    mat::AbstractMatrix{Bool},
)

    S = size(mat, 1)

    pred_has_links = vec(sum(mat, dims=2) .> 0)
    prey_has_links = vec(sum(mat, dims=1) .> 0)

    active = sum(pred_has_links .|| prey_has_links)

    connectance = sum(mat) / (S^2)

    return active, connectance

end


"""
    _rand_categorical(weights)

Sample an index from a probability vector.

The input should already be normalised to sum to one.
"""
function _rand_categorical(weights::AbstractVector{<:Real})

    r = rand()

    cumulative = zero(eltype(weights))

    for i in eachindex(weights)

        cumulative += weights[i]

        if r <= cumulative
            return i
        end

    end

    return lastindex(weights)

end

"""
    _normalise_probability_matrix!(P)

Scale a probability matrix to the interval [0,1].
"""
function _normalise_probability_matrix!(P)

    maxval = maximum(P)

    if maxval > 0 && isfinite(maxval)
        P ./= maxval
    else
        fill!(P, 0.0)
    end

    clamp!(P, 0.0, 1.0)

    return P

end

"""
    _sample_links(matrix, probabilities)

Perform probabilistic link retention using supplied probabilities.

The original matrix is used as a mask so that downsampling can only remove
existing interactions and never create new ones.
"""
function _sample_links(
    matrix::AbstractMatrix{Bool},
    probabilities::AbstractMatrix{<:Real},
)

    return matrix .& (rand(size(matrix)) .<= probabilities)

end

"""
    _existing_links(matrix)

Return the CartesianIndices corresponding to existing interactions.
"""
_existing_links(matrix::AbstractMatrix{Bool}) = findall(matrix)