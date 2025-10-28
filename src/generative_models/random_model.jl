"""
random_model.jl
-----------------
Generates a random food web (Erdős-Rényi model) using Graphs.jl.

This file contains:
1.  `generate_random_model`: A standardised function called by `main.jl`.
    It now filters the result based on the emergent basal species percentage.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# using Graphs, Random

"""
    generate_random_model(S::Int, C::Float64)

Generates a random directed (Erdős-Rényi) graph with `S` nodes
and `L = C * S^2` edges.
Filters the generated web against the global `BASAL_RANGE`.

# Arguments
- `S`: Number of species.
- `C`: Target connectance (L/S^2). This is *randomly drawn* in `main.jl`.

# Returns
- A `NamedTuple` with `adj` and `percent_basal` or `nothing` if filter fails.
"""
function generate_random_model(S::Int, C::Float64)
    
    # 1. Calculate target number of links
    L = round(Int, C * S^2)
    
    # 2. Create the graph
    g = SimpleDiGraph(S, L)
    
    # 3. Get the adjacency matrix
    # adj[i, j] = 1 means edge i -> j (i is predator)
    adj = Matrix(adjacency_matrix(g))
    
    # 4. Calculate emergent % basal species
    # (An outdegree of 0 is a producer)
    producers = findall(i -> outdegree(g, i) == 0, 1:S)
    percent_basal = length(producers) / S
    
    # 5. New Requirement: Filter by emergent basal range
    if !is_in_basal_range(percent_basal)
        # @warn "Random web failed basal range check ($percent_basal). Skipping."
        return nothing
    end

    return (adj=adj, percent_basal=percent_basal)
end