using LinearAlgebra

"""
    trophic_level(A::AbstractMatrix{Bool}; species=nothing, exclude_cannibalism=true)

Calculate species trophic levels using the Williams & Martinez (2004)
shortest-path/diet matrix formulation.

The adjacency matrix convention is:

`A[i,j] == true` means species `i` consumes species `j`.

Basal species (species without prey) receive trophic level 1.

If `exclude_cannibalism=true`, self-consumption links (`A[i,i] == true`)
are removed before calculating trophic levels. This is done by setting the
diagonal of a copy of the adjacency matrix to `false`, leaving the original
input matrix unchanged.

Species that only have cannibalistic links will therefore be treated as
having no prey and will receive a trophic level of 1.

If `species` is provided, returns a dictionary mapping species identifiers
to trophic levels.
"""
function trophic_level(
    A::AbstractMatrix{Bool};
    species=nothing,
    exclude_cannibalism=true
)

    A = copy(A)

    if exclude_cannibalism
        A[diagind(A)] .= false
    end

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

    isnothing(species) ? tls : Dict(zip(species,tls))
end