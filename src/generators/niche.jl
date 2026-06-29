#=
niche_model.jl
Generates a food web using the Niche Model (Williams & Martinez, 2000).

*** MODIFIED VERSION ***
This script is a "pure generator". It does NOT perform any internal
filtering for basal percentage or connectance. All filtering is
handled externally by the main script (01_BuildNetworks.jl).
=#

# --- 1. Dependencies ---
using Distributions, Random 

"""
    niche_model(S::Int, C_target::Float64)

Generates an adjacency matrix using the standard Niche Model algorithm.

*** NOTE: This function does NOT filter the output. ***

# Arguments
- `S::Int`: Species richness.
- `C_target::Float64`: Target connectance, used to parameterize
  the feeding range distribution.

# Returns
- `Matrix{Int}`: The (S x S) adjacency matrix.
- `nothing`: If parameter calculation fails (e.g., C_target > 0.5).
"""
function niche_model(S::Int, C_target::Float64)
    # --- Initialize Adjacency Matrix ---
    adj = zeros(Int, S, S)
    
    # --- 1. Assign Niche Values (n) ---
    n = rand(Uniform(0, 1), S)
    
    # --- 2. Assign Feeding Ranges (r) ---
    beta_param = 1.0 / (2.0 * C_target) - 1.0
    
    # If beta_param is invalid, this model fails. Return `nothing`
    # so the external wrapper function knows it failed.
    if beta_param <= 0
        return nothing
    end
    r = n .* rand(Beta(1.0, beta_param), S)
    
    # --- 3. Assign Feeding Centres (c) ---
    c = zeros(Float64, S)
    for i = 1:S
        lower_c_bound = r[i] / 2.0
        upper_c_bound = n[i] - r[i] / 2.0
        if lower_c_bound >= upper_c_bound
             c[i] = rand(Uniform(0, n[i]))
        else
            c[i] = rand(Uniform(lower_c_bound, upper_c_bound))
        end
    end

    # --- 4. Create Links ---
    for i in 1:S
        lower_bound = c[i] - (r[i] / 2.0)
        upper_bound = c[i] + (r[i] / 2.0)
        for j in 1:S
            if i == j; continue; end
            if lower_bound < n[j] < upper_bound
                adj[i, j] = 1 # Predator i eats Prey j
            end
        end # end prey loop
    end # end predator loop
    
    # --- 5. (DELETED) Filters ---
    # All filtering is now handled by the external wrapper
    # in 01_BuildNetworks.jl
    
    # --- 6. Return Success ---
    # Return the generated matrix directly.
    return adj
end # end generate_niche_model