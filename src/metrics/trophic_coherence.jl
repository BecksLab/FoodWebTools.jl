"""
    trophic_coherence(A::AbstractMatrix{Bool})

Calculate trophic incoherence q.

A[i,j] == true indicates species i consumes species j.

q is the root mean square deviation of trophic distances
from the expected predator-prey separation of one trophic level.

q = 0 indicates perfect trophic coherence.

Returns NaN for networks without interactions.
"""
function trophic_coherence(A::AbstractMatrix{Bool})

    tl = trophic_level(A)

    edges = findall(A)

    isempty(edges) && return NaN

    distances = [
        tl[i] - tl[j]
        for edge in edges
        for (i,j) in [Tuple(edge)]
    ]

    return sqrt(mean((distances .- 1).^2))
end