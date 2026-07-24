# =============================================================================
# Random Downsampling
# =============================================================================

"""
    _single_step(matrix, ::RandomSampling)

Perform a single random downsampling step.

Each existing interaction has an equal probability of being retained.
"""
function _single_step(
    matrix::AbstractMatrix{Bool},
    ::RandomSampling,
)

    # Generate a random retention probability for every existing link.
    #
    # Since the matrix is later masked with the original interaction
    # matrix, probabilities for absent links are irrelevant.
    probabilities = rand(Float64, size(matrix)...)

    return _sample_links(matrix, probabilities)

end


# =============================================================================
# Iterative removal weights
# =============================================================================

"""
    _link_removal_weights(matrix, method::RandomSampling)

Return the candidate links and their removal weights.

For random downsampling every interaction is equally likely to be removed,
so each existing link receives the same weight.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    ::RandomSampling,
)

    links = _existing_links(matrix)

    weights = ones(Float64, length(links))

    return links, weights

end