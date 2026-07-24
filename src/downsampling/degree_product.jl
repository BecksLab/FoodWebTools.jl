# =============================================================================
# Degree-product Downsampling
# =============================================================================

"""
    _degree_product_weights(matrix)

Calculate degree-product scores for existing interactions.

For each interaction:

    consumer → resource

the removal weight is:

    consumer out-degree × resource in-degree

where:

- out-degree = consumer generality (number of resources consumed)
- in-degree  = resource vulnerability (number of consumers)

Interactions involving highly connected species therefore receive higher
removal probability.
"""
function _degree_product_weights(
    matrix::AbstractMatrix{Bool},
)

    links = _existing_links(matrix)

    isempty(links) && return links, Float64[]

    # Consumer generality.
    out_degree = vec(sum(matrix, dims=2))

    # Resource vulnerability.
    in_degree = vec(sum(matrix, dims=1))

    weights = Vector{Float64}(undef, length(links))

    for (idx, link) in enumerate(links)

        consumer = link[1]
        resource = link[2]

        weights[idx] =
            Float64(out_degree[consumer] * in_degree[resource])

    end

    return links, weights

end


# =============================================================================
# Single-step downsampling
# =============================================================================

"""
    _single_step(matrix, method::DegreeProduct)

Perform a single probabilistic degree-product downsampling step.

Links are retained with probability inversely proportional to their
degree-product score.
"""
function _single_step(
    matrix::AbstractMatrix{Bool},
    ::DegreeProduct,
)

    links, weights = _degree_product_weights(matrix)

    isempty(links) && return matrix

    # Convert removal weights into retention probabilities.
    #
    # Large degree-product values correspond to fragile links and therefore
    # receive lower retention probability.
    max_weight = maximum(weights)

    retention = 1 .- (weights ./ max_weight)

    retention .= clamp.(retention, 0.0, 1.0)

    result = copy(matrix)

    for (idx, link) in enumerate(links)

        if rand() > retention[idx]
            result[link] = false
        end

    end

    return result

end


# =============================================================================
# Iterative removal weights
# =============================================================================

"""
    _link_removal_weights(matrix, method::DegreeProduct)

Return removal weights for every existing interaction.

Higher degree-product links are preferentially removed.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    ::DegreeProduct,
)

    return _degree_product_weights(matrix)

end