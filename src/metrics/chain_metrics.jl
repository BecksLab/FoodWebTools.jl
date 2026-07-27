"""
    compute_reachable_to_top(preds, top_set)

Identify all species that can reach at least one top predator.

Given a vector of predecessor lists (`preds`) and the set of top predators,
this function iteratively propagates reachability backwards through the food
web until no additional species can be added.

This is an internal helper used by [`chain_metrics`](@ref).
"""
function compute_reachable_to_top(preds, top_set)

    reachable = Set(top_set)

    changed = true
    while changed
        changed = false

        for s in eachindex(preds)
            if any(in(reachable), preds[s])
                if s ∉ reachable
                    push!(reachable, s)
                    changed = true
                end
            end
        end
    end

    return reachable
end

function weakly_connected_components(A::AbstractMatrix{Bool})

    n = size(A,1)
    visited = falses(n)
    n_components = 0

    for start in 1:n

        # ignore isolated species
        if visited[start] || !(any(@view A[start,:]) || any(@view A[:,start]))
            continue
        end

        n_components += 1

        stack = [start]
        visited[start] = true

        while !isempty(stack)

            u = pop!(stack)

            for v in 1:n

                if !visited[v] && (A[u,v] || A[v,u])
                    visited[v] = true
                    push!(stack, v)
                end

            end
        end
    end

    return n_components
end


"""
    chain_metrics(A::AbstractMatrix{Bool}; max_depth=size(A,1)-1)

Calculate food-chain statistics for a directed food web.

The input is interpreted as a directed adjacency matrix, where

    A[i,j] = true

indicates that species `i` is consumed by species `j` (prey → predator).

A food chain is defined as a simple directed path from a basal species to a top
predator containing at least one trophic interaction. Species may not appear
more than once within a chain, preventing infinite paths in food webs
containing feeding loops or cannibalism.

Basal species are active species with no prey (generality = 0), while top
predators are active species with no predators (vulnerability = 0). Isolated
species (those with no trophic interactions) are excluded from the analysis.

Food chains are enumerated using a depth-first search. The optional
`max_depth` argument limits the maximum chain length explored and may be useful
for very large or highly connected food webs. By default, the search depth is
`size(A,1)-1`, the longest possible simple path.

If the network contains multiple weakly connected components (ignoring edge
direction), food-chain statistics are computed across all components and a
warning is issued.

# Arguments

- `A::AbstractMatrix{Bool}`: directed adjacency matrix.
- `max_depth::Int=size(A,1)-1`: maximum chain length explored.

# Returns

A named tuple containing

- `ChLen::Float64`: mean food-chain length.
- `ChSD::Float64`: standard deviation of food-chain lengths. Returns `0.0`
  when only a single food chain is present.
- `ChNum::Float64`: natural logarithm of the number of food chains.

If the network contains no food chains, `ChLen` and `ChSD` are returned as
`NaN`, while `ChNum` is returned as `0.0`.

# Examples

```julia
A = Bool[
    0 1 0
    0 0 1
    0 0 0
]

chain_metrics(A)

# (ChLen = 2.0, ChSD = 0.0, ChNum = 0.0)
```
"""
function chain_metrics(A::AbstractMatrix{Bool};
                       max_depth::Int=size(A,1)-1)

    n = size(A,1)

    # Number of prey and predators
    row_sum = vec(sum(A, dims=2))
    col_sum = vec(sum(A, dims=1))

    # Species participating in at least one trophic interaction
    active = (row_sum .+ col_sum) .> 0

    # Basal and top species
    basal = findall((row_sum .== 0) .& active)
    top_set = Set(findall((col_sum .== 0) .& active))

    if isempty(basal) || isempty(top_set)
        return (ChLen = NaN,
                ChSD  = NaN,
                ChNum = 0.0)
    end

    # Predecessor lists (prey of each species)
    preds = [findall(@view A[:,j]) for j in 1:n]

    # Warn if multiple weakly connected components exist
    n_components = weakly_connected_components(A)

    if n_components > 1
        @warn "Network contains $n_components disconnected components. Food-chain statistics are computed across all components."
    end

    # Prune species that cannot reach a top predator
    reachable = compute_reachable_to_top(preds, top_set)

    visited = falses(n)
    all_lengths = Int[]

    function dfs(node::Int, depth::Int)

        if depth > max_depth
            return
        end

        if visited[node]
            return
        end

        if node ∉ reachable
            return
        end

        visited[node] = true

        if node in top_set && depth > 0
            push!(all_lengths, depth)
        end

        for nxt in preds[node]
            dfs(nxt, depth + 1)
        end

        visited[node] = false
    end

    for b in basal
        dfs(b, 0)
    end

    if isempty(all_lengths)
        return (ChLen = NaN,
                ChSD  = NaN,
                ChNum = 0.0)
    end

    return (
        ChLen = mean(all_lengths),
        ChSD  = length(all_lengths) == 1 ? 0.0 : std(all_lengths),
        ChNum = log(length(all_lengths))
    )
end