#=
cascade_model.jl
Generates a food web using the Cascade Model (Cohen & Newman, 1985).

*** MODIFIED VERSION ***
This script is a "pure generator". It does NOT perform any internal
filtering for basal percentage or connectance. All filtering is
handled externally by the main script (01_BuildNetworks.jl).
=#

# --- 1. Dependencies ---
using Random

"""
    cascade_model(S::Int, C_target::Float64)

Generates an adjacency matrix using the standard Cascade Model algorithm.
Links (i -> j) only form if `n[i] > n[j]` (cascade constraint) and occur
with probability `p` derived from `C_target`.

*** NOTE: This function does NOT filter the output. ***

# Arguments
- `S::Int`: Species richness.
- `C_target::Float64`: Target connectance, used to calculate link probability `p`.

# Returns
- `Matrix{Int}`: The (S x S) adjacency matrix.
- `nothing`: If probability calculation fails (e.g., p > 1).
"""
function cascade_model(S::Int, C_target::Float64)
    # --- Initialize Adjacency Matrix ---
    adj = zeros(Int, S, S)
    
    # --- 1. Assign Hierarchy Values (niche_values) ---
    niche_values = rand(S)
    
    # --- 2. Calculate Link Probability (p) ---
    if S <= 1
         p = 0.0
    else
         p = (2.0 * C_target * S) / (S - 1.0)
    end
    
    # If p is invalid, this model fails. Return `nothing` so the
    # external wrapper function knows it failed.
    if p > 1.0
        return nothing
    end

    # --- 3. Create Links ---
    for i in 1:S
        for j in 1:S
            if i == j; continue; end
            
            if niche_values[i] > niche_values[j]
                if rand() < p
                    adj[i, j] = 1 # Predator i eats Prey j
                end
            end
        end # end prey loop
    end # end predator loop
    
    # --- 4. (DELETED) Filters ---
    # All filtering is now handled by the external wrapper
    # in 01_BuildNetworks.jl
    
    # --- 5. Return Success ---
    # Return the generated matrix directly.
    return adj
end # end generate_cascade_model