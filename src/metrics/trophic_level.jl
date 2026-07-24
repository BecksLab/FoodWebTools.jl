using LinearAlgebra

"""
    trophic_level(A::AbstractMatrix{Bool}; species=nothing)

Calculate species trophic levels using the Williams & Martinez (2004)
shortest-path/diet matrix formulation.

The adjacency matrix convention is:

`A[i,j] == true` means species `i` consumes species `j`.

Basal species (species without prey) receive trophic level 1.

If `species` is provided, returns a dictionary mapping species identifiers
to trophic levels.
"""
function trophic_level(
    A::AbstractMatrix{Bool};
    species=nothing
)

    S = size(A,1)

    D = zeros(Float64, S, S)

    for i in 1:S
        prey = findall(A[i,:])

        if !isempty(prey)
            D[i, prey] .= -1 / length(prey)
        end
    end

    D[diagind(D)] .+= 1.0

    b = ones(S)

    tls = try
        D \ b
    catch e
        if e isa LinearAlgebra.SingularException
            pinv(D) * b
        else
            rethrow(e)
        end
    end

    if isnothing(species)
        return tls
    end

    return Dict(zip(species,tls))
end