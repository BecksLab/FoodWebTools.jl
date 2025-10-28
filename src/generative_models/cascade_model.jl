"""
cascade_model.jl
-----------------
Generates a food web using the Cascade Model (Cohen & Newman, 1985).

Contains a single function `generate_cascade_model` called by `main.jl`.
Applies external checks for basal percentage and connectance.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# Need Random for uniform draws.
using Random

"""
    generate_cascade_model(S::Int, C_target::Float64)

Generates an adjacency matrix using the standard Cascade Model algorithm.
Links (i -> j) only form if `n[i] > n[j]` (cascade constraint) and occur
with probability `p` derived from `C_target`. Applies filters for emergent
basal percentage and connectance ranges.

# Arguments
- `S::Int`: Species richness.
- `C_target::Float64`: Target connectance (drawn randomly in `main.jl`).
    Used to calculate link probability `p`.

# Returns
- `NamedTuple`: Contains `adj`, `percent_basal`, `connectance` if checks pass.
- `nothing`: If probability calculation fails or checks fail.

---
### Pre-defined Parameters:
* `S::Int`: Species richness.
* `C_target::Float64`: **Target Connectance**. Used to calculate link
    probability `p`.

### Emergent Properties:
* **Adjacency Matrix:** Emerges from random hierarchy (`niche_values`) and
    probabilistic links subject to the cascade constraint.
* **Percent Basal:** Emergent property, checked against `BASAL_RANGE`.
* **Connectance:** Emergent property, checked against `CONNECTANCE_RANGE`.
* **Body Mass:** Not used by this model.
"""
function generate_cascade_model(S::Int, C_target::Float64)
    # --- Initialize Adjacency Matrix ---
    adj = zeros(Int, S, S)
    
    # --- 1. Assign Hierarchy Values (niche_values) ---
    # Each species gets a random value, establishing a hierarchy.
    niche_values = rand(S)
    
    # --- 2. Calculate Link Probability (p) ---
    # The probability `p` of a link forming *given the cascade constraint is met*.
    # Derivation:
    # Expected Links E[L] = C_target * S^2
    # Number of possible links under cascade constraint = S*(S-1)/2
    # E[L] = p * (S*(S-1)/2)
    # p = E[L] / (S*(S-1)/2) = (C_target * S^2) / (S*(S-1)/2) = 2 * C_target * S / (S - 1)
    if S <= 1 # Avoid division by zero if S=1
         p = 0.0
    else
         p = (2.0 * C_target * S) / (S - 1.0)
    end
    # If the calculated probability is > 1 (due to high C_target), it's invalid.
    # Return failure. `main.jl` will try again with a new C_target.
    if p > 1.0
        # @warn "Cascade model C_target ($C_target) too high, p=$p > 1.0"
        return nothing
    end

    # --- 3. Create Links ---
    # Iterate through potential predators.
    for i in 1:S
        # Iterate through potential prey.
        for j in 1:S
            # Prevent cannibalism.
            if i == j; continue; end
            
            # Apply Cascade Constraint: Link i -> j only possible if n[i] > n[j].
            if niche_values[i] > niche_values[j]
                # If constraint met, add link with probability p.
                if rand() < p
                    adj[i, j] = 1 # Predator i eats Prey j
                end
            end
        end # end prey loop
    end # end predator loop
    
    # --- 4. Apply Filters ---

    # --- CHECK 1: Emergent Basal % ---
    percent_basal = calculate_emergent_producers(adj)
    if !is_in_basal_range(percent_basal)
        return nothing # Fail if outside range.
    end

    # --- CHECK 2: Emergent Connectance ---
    connectance = calculate_connectance(adj)
    if !is_in_connectance_range(connectance)
        return nothing # Fail if outside range.
    end
    
    # --- Return Success ---
    # If both checks passed, return the results.
    return (adj=adj, percent_basal=percent_basal, connectance=connectance)
end # end generate_cascade_model