"""
trophic_level(A::AbstractMatrix{Bool};
              species=nothing)

    Calculates the trophic level of all species in a network using the average 
    shortest path from the prey of species 𝑖 to a basal species

    Williams, Richard J., and Neo D. Martinez. 2004. “Limits to Trophic Levels 
    and Omnivory in Complex Food Webs: Theory and Data.” The American Naturalist 
    163 (3): 458–68. https://doi.org/10.1086/381964.
"""
function trophic_level(
    A::AbstractMatrix{Bool};
    species=nothing
)
    # Species richness
    S = size(A, 1)
    in_degree = sum(A; dims = 2)
    # Diet matrix
    D = -(A ./ in_degree)
    D[isnan.(D)] .= 0.0
    D[diagind(D)] .= 1.0 .- D[diagind(D)]
    # Solve with the inverse matrix.
    inverse = iszero(det(D)) ? pinv : inv
    tls = inverse(D) * ones(S)
    tls = vec(inverse(D) * ones(S))

    if isnothing(species)
        return tls
    end

    return Dict(zip(species, tls))
end