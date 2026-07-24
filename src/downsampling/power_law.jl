# =============================================================================
# Power-law Downsampling (Roopnarine, 2006)
# =============================================================================

"""
    _retention_probabilities(matrix, method::PowerLaw)

Compute the relative retention probability for each consumer based on the
power-law scaling model of Roopnarine (2006).

Consumers with larger trophic breadth (higher generality) receive larger
retention probabilities and are therefore less likely to lose interactions.
"""
function _retention_probabilities(
    matrix::AbstractMatrix{Bool},
    method::PowerLaw,
)

    S = size(matrix, 1)

    # Consumer generality (number of prey per consumer).
    generality = vec(sum(matrix, dims=2))

    # Scaling parameter from Roopnarine (2006).
    E = exp(log(S) * (method.y - 1) / method.y)

    probabilities = exp.(generality ./ E)

    # Scale into the interval [0,1].
    probabilities ./= maximum(probabilities)

    return probabilities

end


# =============================================================================
# Single-step downsampling
# =============================================================================

"""
    _single_step(matrix, method::PowerLaw)

Perform a single probabilistic downsampling step using the power-law
retention model.

Every interaction belonging to the same consumer shares the same retention
probability.
"""
function _single_step(
    matrix::AbstractMatrix{Bool},
    method::PowerLaw,
)

    S = size(matrix, 1)

    consumer_probabilities = _retention_probabilities(matrix, method)

    probability_matrix = zeros(Float64, S, S)

    for consumer in 1:S

        prey = findall(matrix[consumer, :])

        probability_matrix[consumer, prey] .= consumer_probabilities[consumer]

    end

    return _sample_links(matrix, probability_matrix)

end


# =============================================================================
# Iterative removal weights
# =============================================================================

"""
    _link_removal_weights(matrix, method::PowerLaw)

Return removal weights for every existing interaction.

Interactions belonging to consumers with low retention probabilities receive
larger removal weights and are therefore more likely to be removed.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    method::PowerLaw,
)

    links = _existing_links(matrix)

    isempty(links) && return links, Float64[]

    consumer_probabilities = _retention_probabilities(matrix, method)

    weights = similar(consumer_probabilities, Float64, length(links))

    for (k, link) in enumerate(links)

        consumer = link[1]

        p_retain = consumer_probabilities[consumer]

        weights[k] = 1.0 - p_retain + 1e-6

    end

    return links, weights

end