"""
random_model.jl
-----------------
Generates a random food web (Erdős-Rényi model variant) using Graphs.jl.

Contains a single function `generate_random_model` called by `main.jl`.
Generates a directed graph with exactly `L = C_target * S^2` links placed randomly.
Applies external checks for basal percentage and connectance.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# Need Graphs for SimpleDiGraph, Random for link placement (internal).
using Graphs, Random

"""
    generate_random_model(S::Int, C_target::Float64)

Generates a random directed graph (Erdős-Rényi (n, L) model) with a fixed
number of links `L` derived from `C_target`. Applies filters for emergent basal
percentage and connectance ranges.

# Arguments
- `S::Int`: Species richness (number of nodes).
- `C_target::Float64`: Target connectance (drawn randomly in `main.jl`).
    Used to calculate the exact number of links `L`.

# Returns
- `NamedTuple`: Contains `adj`, `percent_basal`, `connectance` if checks pass.
- `nothing`: If the connectance check fails (usually due to rounding `L`).

---
### Pre-defined Parameters:
* `S::Int`: Species richness.
* `C_target::Float64`: **Target Connectance**. Used to calculate the exact
    number of links `L = round(C_target * S^2)`.

### Emergent Properties:
* **Adjacency Matrix:** Emerges from the `L` randomly placed directed links.
* **Percent Basal:** Emergent property, checked against `BASAL_RANGE`.
* **Connectance:** Calculated exactly as `L / S^2`. Checked against
    `CONNECTANCE_RANGE`. (Might fail if rounding `L` pushes C just outside).
* **Body Mass:** Not used by this model.
"""
function generate_random_model(S::Int, C_target::Float64)
    
    # --- 1. Calculate Target Number of Links (L) ---
    # Round the target number of links based on C_target and S.
    L = round(Int, C_target * S^2)
    # Ensure L is within valid bounds [0, S*(S-1)].
    L = max(0, min(L, S * (S - 1))) # Max possible links excluding self-loops.

    # --- 2. Create Random Directed Graph ---
    # Use Graphs.jl to create a SimpleDiGraph with S nodes and exactly L edges
    # placed uniformly at random among possible directed edges (excluding self-loops).
    g = SimpleDiGraph(S, L)
    
    # --- 3. Get Adjacency Matrix ---
    # Convert the graph object to its adjacency matrix representation.
    # adj[i, j] = 1 means edge i -> j (Predator i -> Prey j).
    adj = Matrix(adjacency_matrix(g))
    
    # --- 4. Apply Filters ---

    # --- CHECK 1: Emergent Basal % ---
    # Find nodes with out-degree 0 in the graph `g`.
    producers = findall(i -> outdegree(g, i) == 0, 1:S)
    percent_basal = length(producers) / S
    # Check if the calculated percentage is outside the allowed range.
    if !is_in_basal_range(percent_basal)
        return nothing # Fail if outside range.
    end

    # --- CHECK 2: Emergent Connectance ---
    # Calculate connectance exactly based on the fixed L.
    connectance = L / (S^2)
    # Check if this connectance falls within the allowed range.
    # Note: This check might fail if `round(Int, C_target * S^2)` produces an L
    # such that L/S^2 is slightly outside the original C_target range.
    if !is_in_connectance_range(connectance)
        return nothing # Fail if outside range.
    end

    # --- Return Success ---
    # If both checks passed, return the results.
    return (adj=adj, percent_basal=percent_basal, connectance=connectance)
end # end generate_random_model