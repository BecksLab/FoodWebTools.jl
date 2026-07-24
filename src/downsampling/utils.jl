# =============================================================================
# Shared Utility Functions
# =============================================================================

"""
    _active_species(matrix)

Return the number of active species in a binary interaction network.

A species is considered *active* if it participates in at least one realised
interaction, either as a consumer (row) or as a resource (column).

This definition treats species with at least one incoming *or* outgoing
interaction as present in the realised food web.

# Returns

An integer giving the number of active species.
"""
function _active_species(
    matrix::AbstractMatrix{Bool},
)

    consumers = vec(sum(matrix, dims=2) .> 0)
    resources = vec(sum(matrix, dims=1) .> 0)

    return sum(consumers .|| resources)

end


"""
    _connectance(matrix)

Calculate the connectance of a binary interaction network.

Connectance is defined as

    L / S²

where

- `L` is the number of realised interactions, and
- `S` is the number of species.

# Returns

A value between 0 and 1.
"""
function _connectance(
    matrix::AbstractMatrix{Bool},
)

    S = size(matrix, 1)

    return sum(matrix) / (S^2)

end


"""
    _rand_categorical(weights)

Sample an index from a probability vector using inverse-CDF sampling.

The supplied weights should already be normalised so that they sum to one.

# Returns

The index of the sampled element.
"""
function _rand_categorical(
    weights::AbstractVector{<:Real},
)

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

Normalise a matrix of probabilities to the interval `[0,1]`.

The matrix is scaled by its maximum finite value before being clamped to the
valid probability range. If all entries are zero (or non-finite), the matrix
is filled with zeros.

The operation is performed in-place.
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

Perform a single probabilistic downsampling step.

Each realised interaction is retained independently according to its
corresponding probability in `probabilities`.

The original interaction matrix is used as a mask so that downsampling can
remove existing links but can never create new interactions.
"""
function _sample_links(
    matrix::AbstractMatrix{Bool},
    probabilities::AbstractMatrix{<:Real},
)

    draws = rand(Float64, size(matrix))

    return matrix .& (draws .<= probabilities)

end


"""
    _existing_links(matrix)

Return the locations of all realised interactions.

# Returns

A vector of `CartesianIndex{2}` objects corresponding to entries equal to
`true` in the interaction matrix.
"""
_existing_links(matrix::AbstractMatrix{Bool}) = findall(matrix)


# =============================================================================
# Iterative downsampling helpers
# =============================================================================

"""
    _valid_removal(matrix, link; allow_species_loss, min_species)

Determine whether removing a candidate interaction satisfies the ecological
constraints imposed by the iterative downsampling algorithm.

The candidate removal is accepted if

1. the resulting network retains at least `min_species` active species, and
2. if `allow_species_loss == false`, no species becomes disconnected as a
   direct consequence of the removal.

Returns `true` if the removal is valid and `false` otherwise.
"""
function _valid_removal(
    matrix::AbstractMatrix{Bool},
    link::CartesianIndex;
    allow_species_loss::Bool,
    min_species::Int,
)

    temp = copy(matrix)
    temp[link] = false

    active = _active_species(temp)

    active < min_species && return false

    if !allow_species_loss

        before = _active_species(matrix)

        if active < before
            return false
        end

    end

    return true

end


"""
    _choose_removal(matrix, links, weights; kwargs...)

Sample a valid interaction for removal.

Candidate links are sampled according to the supplied removal weights. If a
sampled link violates the ecological constraints enforced by
`_valid_removal`, it is discarded and another candidate is sampled from the
remaining links.

This continues until either

- a valid link is found, or
- no valid candidates remain.

# Returns

A `CartesianIndex` identifying the chosen interaction, or `nothing` if no
valid removal exists.
"""
function _choose_removal(
    matrix::AbstractMatrix{Bool},
    links,
    weights;
    allow_species_loss::Bool,
    min_species::Int,
)

    candidates = collect(links)
    probs = Float64.(weights)

    while !isempty(candidates)

        probs ./= sum(probs)

        chosen = _rand_categorical(probs)

        link = candidates[chosen]

        if _valid_removal(
            matrix,
            link;
            allow_species_loss,
            min_species,
        )
            return link
        end

        deleteat!(candidates, chosen)
        deleteat!(probs, chosen)

    end

    return nothing

end