"""
ltm.jl
-----------------
Generates a food web using the Latent Trait Model (LTM).

Contains the core `ltm` generation function (provided by user) and a
standardised wrapper function `generate_ltm_model` called by `main.jl`.

This version runs the LTM deterministically and without internal viability checks,
as requested. External checks for basal percentage and connectance are applied
in the wrapper.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# Load necessary packages for the LTM calculations.
using Distributions, Graphs, Random, Statistics

"""
# --- 2. LTM Generator (User Provided Code) ---
    ltm(...) -> NamedTuple

Generates a food web interaction matrix using a Latent Trait Model (LTM). This
version integrates correlated latent trait generation and offers highly
flexible control over structural viability checks.

The function models the probability of a feeding link between a predator and a
prey based on two main components: their body mass ratio and their positions
along two latent trait axes (vulnerability and foraging).

### Core Features:
1.  **Correlated Trait Generation:** Instead of requiring pre-generated traits, this
    function creates them internally using a multivariate normal distribution. This
    allows the user to specify a precise Pearson correlation coefficient (`ρ`)
    between the foraging and vulnerability traits, enabling the simulation of
    ecological trade-offs (e.g., "live fast, die young").
2.  **Flexible Viability Workflow:** The structural viability of the generated
    network can be handled in three distinct ways, controlled by the `ensure_viability`
    and `remove_flawed` flags, making the function suitable for a wide range of
    exploratory and validation tasks.
3.  **Stochastic & Deterministic Modes:** Can generate webs probabilistically,
    allowing for variation, or deterministically, creating a web with a fixed
    number of links based on the highest interaction probabilities.

# Arguments
- `species_indices::AbstractVector{Int}`: A vector of integer indices for the
  species to be included (e.g., `1:S`).
- `bodymasses::AbstractVector{Float64}`: A vector of body masses for each species.
  Crucial for calculating the body mass ratio component of the model.
- `metabolic_classes::AbstractVector{Symbol}`: A vector of metabolic classes
  (e.g., `:producer`, `:consumer`). Species marked as `:producer` are prevented
  from having any foraging links.

# Keyword Arguments
## ltm Parameters
- `P_max_target::Float64 = 0.85`: The target maximum probability of an interaction
  at the optimal body mass ratio. Influences the overall density of the web.
- `x_opt::Float64 = 100.0`: The optimal predator-prey body mass ratio (prey mass /
  predator mass) where the interaction probability is highest.
- `sigma_x::Float64 = 1.5`: The standard deviation of the feeding niche along the
  log10 body mass ratio axis. Controls the degree of size specialisation.
- `δ::Float64 = 8.0`: The strength of the latent trait matching. This scales the
  overall influence of the `v_trait * f_trait` interaction term.

## Latent Trait Generation
- `trait_sd::Float64 = 1.0`: The standard deviation for the latent traits.
  Controls the spread of trait values in the community.
- `trait_correlation::Float64 = 0.0`: The desired Pearson correlation coefficient
  (ρ) between vulnerability and foraging traits. Can range from -1.0 to 1.0.

## Generation & Viability Control
- `stochastic::Bool = true`: If `true`, links are formed probabilistically. If
  `false`, links are assigned deterministically to the pairs with the highest
  probabilities until the expected number of links is met.
- `ensure_viability::Bool = false`: The master switch for viability checking.
    - If `false`, the function performs no checks and returns the first web it
      generates with a status of `:unchecked`.
    - If `true`, the function will check the web against the `is_viable` rules.
- `remove_flawed::Bool = false`: This flag is **only** active when `ensure_viability = true`.
    - If `false`, and the generated web is not viable, the function will immediately
      return the flawed web with a status of `:structurally_flawed`.
    - If `true`, and the web is not viable, the function will discard it and
      attempt to generate a new one, up to `max_attempts_internal` times.

# Returns
- A `NamedTuple` containing the following fields:
    - `status::Symbol`: A symbol indicating the outcome. Can be one of:
        - `:unchecked`: `ensure_viability` was `false`.
        - `:viable`: `ensure_viability` was `true` and the web passed the check.
        - `:structurally_flawed`: `ensure_viability` was `true`, `remove_flawed` was `false`, and the web failed the check.
        - `:generation_failed`: `ensure_viability` and `remove_flawed` were `true`, but no viable web could be found in the allowed attempts.
    - `binary_matrix::Matrix{Int}`: The final (S x S) adjacency matrix of the food web.
    - `probability_matrix::Matrix{Float64}`: The (S x S) matrix of calculated interaction probabilities.
    - `v_traits::Vector{Float64}`: The generated vulnerability traits for each species.
    - `f_traits::Vector{Float64}`: The generated foraging traits for each species.
This flexible version can EITHER:
1. Generate correlated traits internally given `trait_sd` and `trait_correlation`.
2. Accept pre-generated `vulnerability_traits` and `foraging_traits`.
"""
function ltm(
    species_indices::AbstractVector{Int},
    bodymasses::AbstractVector{Float64},
    metabolic_classes::AbstractVector{Symbol}; # Expects :producer or :invertebrate
    # LTM Parameters
    P_max_target::Float64 = 0.50,
    x_opt::Float64 = 100.0, # Optimal predator-prey mass ratio (Ropt/Z)
    sigma_x::Float64 = 1.5,
    δ::Float64 = 3.0,
    # Trait Generation
    trait_sd::Float64 = 1.0,
    trait_correlation::Float64 = 0.0,
    vulnerability_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    foraging_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    # Generation Settings (Updated Defaults)
    stochastic::Bool = false,        # Deterministic mode
    ensure_viability::Bool = false,  # Skip internal viability check
    remove_flawed::Bool = false
)
    # --- Initialisation ---
    S = length(species_indices)
    local f_traits, v_traits # Define in local scope for conditional assignment

    # --- BLOCK 1: Generate or Use Provided Latent Traits ---
    # Check if traits need to be generated internally.
    if isnothing(vulnerability_traits) || isnothing(foraging_traits)
        # --- Internal Trait Generation ---
        # Define mean vector (centered at 0).
        μ = [0.0, 0.0]
        # Define covariance matrix based on SD and correlation.
        Σ = [trait_sd^2               trait_correlation*trait_sd^2;
             trait_correlation*trait_sd^2 trait_sd^2]
        # Create multivariate normal distribution.
        dist = MvNormal(μ, Σ)
        # Sample traits (2xS matrix).
        traits_matrix = rand(dist, S)
        f_traits = traits_matrix[1, :] # Row 1: Foraging
        v_traits = traits_matrix[2, :] # Row 2: Vulnerability
        
        # --- Apply Producer Rule ---
        # Find indices of producers.
        producer_indices = findall(mc -> mc == :producer, metabolic_classes)
        # Set foraging trait of producers to 0.
        f_traits[producer_indices] .= 0.0
    else
        # --- Use Provided Traits ---
        v_traits = vulnerability_traits
        f_traits = foraging_traits
        # Assume caller has correctly handled producer foraging traits.
    end # end trait generation block

    # --- BLOCK 2: Derive LTM Coefficients ---
    # Convert ecologically intuitive parameters into quadratic equation coefficients.
    log_x_opt = log10(x_opt)
    _gamma = -1.0 / (2.0 * sigma_x^2)
    _beta = -2.0 * _gamma * log_x_opt
    # Use logit transform for the target max probability.
    logit_P_max_target = log(P_max_target / (1.0 - P_max_target))
    _alpha = logit_P_max_target - (_beta * log_x_opt) - (_gamma * log_x_opt^2)

    # --- BLOCK 3: Generate Interaction Probabilities ---
    # Initialize matrix to store link probabilities.
    prob_matrix = zeros(Float64, S, S)
    # Iterate through all potential predator-prey pairs.
    for pred_idx in species_indices
        # Skip if the potential predator is a producer or has invalid body mass.
        if metabolic_classes[pred_idx] == :producer || bodymasses[pred_idx] <= 0; continue; end
        for prey_idx in species_indices
            # Skip self-loops (cannibalism) or invalid prey body mass.
            if pred_idx == prey_idx || bodymasses[prey_idx] <= 0; continue; end
            
            # Calculate log10 of mass ratio (Prey / Predator).
            log_mass_ratio = log10(bodymasses[prey_idx] / bodymasses[pred_idx])
            # Skip if mass ratio calculation resulted in non-finite value (e.g., log(0)).
            if !isfinite(log_mass_ratio); continue; end
            
            # Calculate body-size component of interaction probability (log-odds).
            body_size_log_odds = _alpha + (_beta * log_mass_ratio) + (_gamma * log_mass_ratio^2)
            # Calculate latent-trait component (log-odds).
            latent_trait_log_odds = δ * v_traits[prey_idx] * f_traits[pred_idx]
            # Combine components.
            total_log_odds = body_size_log_odds + latent_trait_log_odds
            
            # Convert total log-odds back to probability using the logistic function.
            prob_matrix[pred_idx, prey_idx] = 1.0 / (1.0 + exp(-total_log_odds))
        end # end prey loop
    end # end predator loop

    # --- BLOCK 4: Build Binary Adjacency Matrix ---
    # Initialize the binary matrix (0s initially).
    binary_matrix = zeros(Int, S, S)
    # Check if stochastic or deterministic mode is selected.
    if stochastic # (Will be false based on updated defaults)
        # --- Stochastic Mode ---
        # For each possible link, compare a random number to the probability.
        for i in 1:S, j in 1:S
            if rand() < prob_matrix[i, j]; binary_matrix[i, j] = 1; end
        end
    else
        # --- Deterministic Mode ---
        # Calculate the expected number of links by summing probabilities.
        expected_links = round(Int, sum(filter(isfinite, prob_matrix)))
        # Proceed only if expected links is positive and finite.
        if expected_links > 0 && isfinite(expected_links)
            # Determine how many links to actually add (cannot exceed S*S).
            # Get indices of all finite probabilities in the flattened matrix.
            finite_indices = findall(isfinite, vec(prob_matrix))
            num_finite = length(finite_indices)
            # Ensure we don't try to take more links than available finite probabilities.
            num_to_take = min(expected_links, num_finite)
            
            # Find the indices corresponding to the `num_to_take` highest probabilities.
            if num_to_take > 0
                # Get the finite probabilities themselves.
                finite_probs = vec(prob_matrix)[finite_indices]
                # Find the permutation that sorts these probabilities in descending order.
                p = sortperm(finite_probs, rev=true)
                # Select the indices within the `finite_indices` array corresponding to the top probabilities.
                top_indices_relative = p[1:num_to_take]
                # Map these relative indices back to the original flattened matrix indices.
                top_indices = finite_indices[top_indices_relative]
                # Set the corresponding positions in the binary matrix to 1.
                binary_matrix[top_indices] .= 1
            end # end if num_to_take > 0
        end # end if expected_links > 0
    end # end deterministic mode

    # --- BLOCK 5: Return Result ---
    # Since ensure_viability is now false by default, we skip internal checks.
    # The status is always :unchecked in this configuration.
    return (status=:unchecked, binary_matrix=binary_matrix, probability_matrix=prob_matrix, v_traits=v_traits, f_traits=f_traits)

end # End of ltm function

# --- 3. Viability Checker (User Provided Code - Not used by default now) ---
"""
    is_viable(g::SimpleDiGraph) -> Bool

A unified, robust function to check if a food web graph is structurally viable
based on its emergent properties. This is a default checker that can be supplied
to the ltm function.

A web is viable if:
1. It is not empty.
2. It contains at least one producer (a species with no prey).
3. All producers are consumed by at least one other species.
4. No species are completely isolated (i.e., have no links at all).
"""
function is_viable(g::SimpleDiGraph)
    S = nv(g)
    if S == 0; return false; end
    producers = findall(i -> outdegree(g, i) == 0, 1:S)
    if isempty(producers); return false; end
    for p in producers
        if indegree(g, p) == 0; return false; end
    end
    for i in 1:S
        if degree(g, i) == 0; return false; end
    end
    return true
end # end is_viable

