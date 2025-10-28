"""
niche_model.jl
-----------------
Generates a food web using the Niche Model (Williams & Martinez, 2000).

Contains a single function `generate_niche_model` called by `main.jl`.
Applies external checks for basal percentage and connectance.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# Need Distributions for Beta distribution, Random for uniform draws.
using Distributions, Random 

"""
    generate_niche_model(S::Int, C_target::Float64)

Generates an adjacency matrix using the standard Niche Model algorithm.
Applies filters for emergent basal percentage and connectance ranges.

# Arguments
- `S::Int`: Species richness.
- `C_target::Float64`: Target connectance (drawn randomly in `main.jl`).
    Used to parameterize the feeding range distribution.

# Returns
- `NamedTuple`: Contains `adj`, `percent_basal`, `connectance` if checks pass.
- `nothing`: If parameter calculation fails or checks fail.

---
### Pre-defined Parameters:
* `S::Int`: Species richness.
* `C_target::Float64`: **Target Connectance**. Used to parameterise the `Beta`
    distribution for feeding ranges `r`.

### Emergent Properties:
* **Adjacency Matrix:** Emerges from random niche values (`n`), ranges (`r`),
    and centers (`c`).
* **Percent Basal:** Emergent property, checked against `BASAL_RANGE`.
* **Connectance:** Emergent property, checked against `CONNECTANCE_RANGE`.
* **Body Mass:** Not used by this model.
"""
function generate_niche_model(S::Int, C_target::Float64)
    # --- Initialize Adjacency Matrix ---
    adj = zeros(Int, S, S)
    
    # --- 1. Assign Niche Values (n) ---
    # Each species gets a random niche value between 0 and 1.
    n = rand(Uniform(0, 1), S)
    
    # --- 2. Assign Feeding Ranges (r) ---
    # Calculate the Beta distribution parameter `beta_param` needed to achieve
    # the target connectance `C_target` on average.
    # Formula: beta = 1/(2C) - 1
    beta_param = 1.0 / (2.0 * C_target) - 1.0
    # If C_target is too high (>0.5), beta_param becomes non-positive,
    # which is invalid for the Beta distribution. Return failure.
    if beta_param <= 0
        # @warn "Niche model C_target ($C_target) too high, beta_param=$beta_param <= 0"
        return nothing
    end
    # Draw feeding range `r` for each species: r = n * Beta(1, beta_param)
    # This ensures larger-niche species tend to have larger ranges.
    r = n .* rand(Beta(1.0, beta_param), S)
    
    # --- 3. Assign Feeding Centres (c) ---
    # Initialize centre vector.
    c = zeros(Float64, S)
    # For each species, draw a feeding centre `c`.
    for i = 1:S
        # The centre must be placed such that the range [c - r/2, c + r/2]
        # falls entirely within [0, n]. The valid interval for c is [r/2, n - r/2].
        # However, we must also ensure the upper bound doesn't exceed 1.0 if n is close to 1.
        lower_c_bound = r[i] / 2.0
        upper_c_bound = n[i] - r[i] / 2.0 # Initial upper bound
        # Handle cases where the interval is invalid (r > n) or very small.
        if lower_c_bound >= upper_c_bound
            # If r >= n, the valid interval is [0, n]. Draw c from this.
             c[i] = rand(Uniform(0, n[i]))
        else
            # Draw c from the valid interval [r/2, n - r/2].
            c[i] = rand(Uniform(lower_c_bound, upper_c_bound))
        end
    end

    # --- 4. Create Links ---
    # Iterate through potential predators.
    for i in 1:S
        # Define the lower and upper bounds of predator i's feeding niche.
        lower_bound = c[i] - (r[i] / 2.0)
        upper_bound = c[i] + (r[i] / 2.0)
        # Iterate through potential prey.
        for j in 1:S
            # Prevent cannibalism.
            if i == j; continue; end
            # Add a link if the prey's niche value `n[j]` falls within the predator's range.
            if lower_bound < n[j] < upper_bound
                adj[i, j] = 1 # Predator i eats Prey j
            end
        end # end prey loop
    end # end predator loop
    
    # --- 5. Apply Filters ---
    
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
end # end generate_niche_model