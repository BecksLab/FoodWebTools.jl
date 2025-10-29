#=
random_model.jl
Generates a random food web (Erdős-Rényi model variant).

*** MODIFIED VERSION ***
This script is a "pure generator". It does NOT perform any internal
filtering for basal percentage or connectance. All filtering is
handled externally by the main script (01_BuildNetworks.jl).
=#

# --- 1. Dependencies ---
using Graphs, Random

"""
    generate_random_model(S::Int, C_target::Float64)

Generates a random directed graph (Erdős-Rényi (n, L) model) with a fixed
number of links `L` derived from `C_target`.

*** NOTE: This function does NOT filter the output. ***

# Arguments
- `S::Int`: Species richness (number of nodes).
- `C_target::Float64`: Target connectance, used to calculate `L`.

# Returns
- `Matrix{Int}`: The (S x S) adjacency matrix.
"""
function generate_random_model(S::Int, C_target::Float64)
    
    # --- 1. Calculate Target Number of Links (L) ---
    L = round(Int, C_target * S^2)
    # Ensure L is within valid bounds [0, S*(S-1)].
    L = max(0, min(L, S * (S - 1))) # Max possible links excluding self-loops.

    # --- 2. Create Random Directed Graph ---
    g = SimpleDiGraph(S, L)
    
    # --- 3. Get Adjacency Matrix ---
    adj = Matrix(adjacency_matrix(g))
    
    # --- 4. (DELETED) Filters ---
    # All filtering is now handled by the external wrapper
    # in 01_BuildNetworks.jl

    # --- 5. Return Success ---
    # Return the generated matrix directly.
    return adj
end # end generate_random_model