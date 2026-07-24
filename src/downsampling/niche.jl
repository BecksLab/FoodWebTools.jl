# =============================================================================
# Topological Niche Downsampling
# =============================================================================

"""
    _get_svd_niche_positions(matrix)

Estimate a one-dimensional niche axis from the interaction matrix using the
leading singular vectors of the Singular Value Decomposition (SVD).

Species are assigned positions along a common latent niche dimension, which
is subsequently used to estimate consumer niche centres and niche breadths.
"""
function _get_svd_niche_positions(matrix::AbstractMatrix{Bool})

    S = size(matrix, 1)

    svd_decomp = svd(Float64.(matrix))

    u1 = svd_decomp.U[:, 1]
    v1 = svd_decomp.V[:, 1]
    s1 = svd_decomp.S[1]

    # Scale coordinates by the dominant singular value.
    consumer_pos = u1 .* sqrt(s1)
    resource_pos = v1 .* sqrt(s1)

    niche = (consumer_pos .+ resource_pos) ./ 2

    lo = minimum(niche)
    hi = maximum(niche)

    if hi > lo
        return (niche .- lo) ./ (hi - lo)
    end

    # Degenerate case: assign evenly spaced positions.
    return collect(range(0.0, 1.0; length=S))

end


# =============================================================================
# Single-step downsampling
# =============================================================================

"""
    _single_step(matrix, method::Niche)

Perform a single probabilistic niche-based downsampling step.

Links closer to the consumer's niche centroid have a higher probability of
being retained.
"""
function _single_step(
    matrix::AbstractMatrix{Bool},
    method::Niche,
)

    S = size(matrix, 1)

    niche_axis = _get_svd_niche_positions(matrix)

    probabilities = zeros(Float64, S, S)

    for consumer in 1:S

        prey = findall(matrix[consumer, :])

        isempty(prey) && continue

        centroid = mean(niche_axis[prey])

        prey_positions = niche_axis[prey]

        σ = length(prey_positions) > 1 ? std(prey_positions) : 0.1
        σ = max(σ, 0.1)

        σ *= method.sigma_scale

        for resource in prey

            distance = abs(niche_axis[resource] - centroid)

            probabilities[consumer, resource] =
                exp(-(distance / σ)^2)

        end

    end

    _normalise_probability_matrix!(probabilities)

    return _sample_links(matrix, probabilities)

end


# =============================================================================
# Iterative removal weights
# =============================================================================

"""
    _link_removal_weights(matrix, method::Niche)

Return removal weights for every existing interaction.

Links lying further from the consumer's niche centre receive larger removal
weights and are therefore more likely to be removed.
"""
function _link_removal_weights(
    matrix::AbstractMatrix{Bool},
    method::Niche,
)

    links = _existing_links(matrix)

    isempty(links) && return links, Float64[]

    niche_axis = _get_svd_niche_positions(matrix)

    weights = Vector{Float64}(undef, length(links))

    for (idx, link) in enumerate(links)

        consumer = link[1]
        resource = link[2]

        prey = findall(matrix[consumer, :])

        centroid = mean(niche_axis[prey])

        prey_positions = niche_axis[prey]

        σ = length(prey_positions) > 1 ? std(prey_positions) : 0.1
        σ = max(σ, 0.1)

        σ *= method.sigma_scale

        distance = abs(niche_axis[resource] - centroid)

        p_retain = exp(-(distance / σ)^2)

        weights[idx] = 1.0 - p_retain + 1e-6

    end

    return links, weights

end