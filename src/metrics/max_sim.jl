"""
    max_sim(A::AbstractMatrix{Bool}) -> Float64

Calculate the mean maximum trophic similarity of a food web.

The input is interpreted as a directed adjacency matrix, where

    A[i,j] = true

indicates that species `i` is consumed by species `j` (prey → predator).

The trophic similarity between two species is defined as

    S(i,j) = (shared prey + shared predators) /
             (total prey + total predators)

where

- shared prey is the number of prey species common to both species,
- shared predators is the number of predators common to both species,
- total prey is the size of the union of their prey,
- total predators is the size of the union of their predators.

For each species, the maximum trophic similarity to every other species is
determined. The reported value is the mean of these maxima across all species.

Species with no trophic interactions in common with any other species
contribute a similarity of `0`.

# Arguments

- `A::AbstractMatrix{Bool}`: directed adjacency matrix.

# Returns

The mean maximum trophic similarity.

# Examples

```julia
A = Bool[
    0 1 1
    0 0 0
    0 0 0
]

max_sim(A)
```
"""
function max_sim(A::AbstractMatrix{Bool})

    A = copy(A)
    A[diagind(A)] .= false

    n = size(A, 1)

    max_similarity = zeros(Float64, n)

    for i in 1:n

        prey_i = @view A[:, i]
        pred_i = @view A[i, :]

        best = 0.0

        for j in 1:n

            i == j && continue

            prey_j = @view A[:, j]
            pred_j = @view A[j, :]

            shared_prey = count(prey_i .& prey_j)
            shared_pred = count(pred_i .& pred_j)

            total_prey = count(prey_i .| prey_j)
            total_pred = count(pred_i .| pred_j)

            denom = total_prey + total_pred

            similarity = denom == 0 ? 0.0 :
                         (shared_prey + shared_pred) / denom

            best = max(best, similarity)

        end

        max_similarity[i] = best

    end

    return mean(max_similarity)

end