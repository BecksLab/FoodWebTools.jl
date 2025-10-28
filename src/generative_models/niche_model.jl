"""
niche_model.jl
-----------------
Generates a food web using the Niche Model (Williams & Martinez, 2000).

This file contains:
1.  `generate_niche_model`: A standardised function called by `main.jl`.
    It now filters the result based on the emergent basal species percentage.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# using Distributions, Random

"""
    generate_niche_model(S::Int, C::Float64)

Generates an adjacency matrix using the standard Niche Model algorithm.
Filters the generated web against the global `BASAL_RANGE`.

# Arguments
- `S`: Number of species.
- `C`: Target connectance (L/S^2). This is *randomly drawn* in `main.jl`.

# Returns
- A `NamedTuple` with `adj` and `percent_basal` or `nothing` if filter fails.
"""
function generate_niche_model(S::Int, C::Float64)
    adj = zeros(Int, S, S)
    
    # 1. Assign niche values `n`
    n = rand(Uniform(0, 1), S)
    
    # 2. Assign feeding ranges `r`
    beta_param = 1.0 / (2.0 * C) - 1.0
    if beta_param <= 0
        # This can happen if C is very high.
        # We return nothing, and main.jl will try again with a new C.
        return nothing
    end
    r = n .* rand(Beta(1.0, beta_param), S)
    
    # 3. Assign feeding centres `c`
    c = zeros(Float64, S)
    for i = 1:S
        # Ensure centre is valid
        c[i] = rand(Uniform(r[i] / 2.0, min(n[i], 1.0 - r[i] / 2.0)))
        # Handle edge case where r > n (rare, but possible)
        if r[i] > n[i]
             c[i] = rand(Uniform(0, n[i]))
        end
    end

    # 4. Create links
    for i in 1:S  # Predator
        lower_bound = c[i] - (r[i] / 2.0)
        upper_bound = c[i] + (r[i] / 2.0)
        
        for j in 1:S  # Prey
            if i == j; continue; end # No cannibalism
            
            if lower_bound < n[j] < upper_bound
                adj[i, j] = 1
            end
        end
    end
    
    # 5. Calculate emergent % basal species
    percent_basal = calculate_emergent_producers(adj)
    
    # 6. New Requirement: Filter by emergent basal range
    if !is_in_basal_range(percent_basal)
        # @warn "Niche web failed basal range check ($percent_basal). Skipping."
        return nothing
    end

    return (adj=adj, percent_basal=percent_basal)
end