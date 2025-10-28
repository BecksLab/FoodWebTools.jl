"""
cascade_model.jl
-----------------
Generates a food web using the Cascade Model (Cohen & Newman, 1985).

This file contains:
1.  `generate_cascade_model`: A standardised function called by `main.jl`.
    It now filters the result based on the emergent basal species percentage.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# using Random

"""
    generate_cascade_model(S::Int, C::Float64)

Generates an adjacency matrix using the standard Cascade Model algorithm.
Filters the generated web against the global `BASAL_RANGE`.

# Arguments
- `S`: Number of species.
- `C`: Target connectance (L/S^2). This is *randomly drawn* in `main.jl`.

# Returns
- A `NamedTuple` with `adj` and `percent_basal` or `nothing` if filter fails.
"""
function generate_cascade_model(S::Int, C::Float64)
    adj = zeros(Int, S, S)
    
    # 1. Assign niche values `n` (species hierarchy)
    niche_values = rand(S)
    
    # 2. Calculate link probability
    # p = (C * S^2) / (S*(S-1)/2) = 2*C*S / (S-1)
    p = (2.0 * C * S) / (S - 1.0)
    if p > 1.0
        # Connectance is too high. Return nothing, main.jl will try again.
        return nothing
    end

    # 3. Create links
    for i in 1:S  # Predator
        for j in 1:S  # Prey
            if i == j; continue; end # No cannibalism
            
            # Rule 1: Predator must be "higher" than prey
            if niche_values[i] > niche_values[j]
                # Rule 2: Add link with probability p
                if rand() < p
                    adj[i, j] = 1
                end
            end
        end
    end
    
    # 4. Calculate emergent % basal species
    percent_basal = calculate_emergent_producers(adj)
    
    # 5. New Requirement: Filter by emergent basal range
    if !is_in_basal_range(percent_basal)
        # @warn "Cascade web failed basal range check ($percent_basal). Skipping."
        return nothing
    end

    return (adj=adj, percent_basal=percent_basal)
end