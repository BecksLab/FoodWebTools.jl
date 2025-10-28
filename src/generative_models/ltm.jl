"""
ltm_model.jl
-----------------
Generates a food web using the Latent Trait Model (LTM).

This file contains:
1.  `is_viable`: The viability checker function.
2.  `LTM`: The main LTM generation function you provided.
3.  `generate_ltm_model`: A standardised wrapper function called by `main.jl`.
    This wrapper now also filters the result based on the emergent
    basal species percentage.
"""

# --- 1. Dependencies (loaded by main.jl) ---
using Distributions, Graphs, Random, Statistics

"""
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
    metabolic_classes::AbstractVector{Symbol};
    # LTM Parameters
    P_max_target::Float64 = 0.85,
    x_opt::Float64 = 100.0,
    sigma_x::Float64 = 1.5,
    δ::Float64 = 8.0,
    # --- MODIFICATION 1: Trait parameters ---
    # These are now optional and will be used for internal generation
    trait_sd::Float64 = 1.0,
    trait_correlation::Float64 = 0.0,
    # Allow for traits to be passed in
    vulnerability_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    foraging_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    # Generation & Viability Control
    stochastic::Bool = true,
    ensure_viability::Bool = false,
    remove_flawed::Bool = false
)
    # Get the number of species (S) from the input indices.
    S = length(species_indices)
    
    # Define trait variables in the local scope so they can be assigned
    # by the conditional logic block below and then used later.
    local f_traits, v_traits

    # --- BLOCK 1: Generate OR Use Provided Latent Traits ---
    # This block is the core fix. It checks if traits were passed in.
    # If not, it generates them. If they were, it uses them.
    # --- MODIFICATION 2: Add conditional logic ---
    if isnothing(vulnerability_traits) || isnothing(foraging_traits)
        # CASE A: Traits were NOT provided. Generate them internally.
        # This branch is used by `setup_recipient_community`.
        # @info "LTM_matrix: No traits provided, generating internally..."

        # 1a. Define the mean vector (μ) for the two traits, centred at 0.0.
        μ = [0.0, 0.0]

        # 1b. Construct the covariance matrix (Σ).
        # Diagonal elements are the variance (σ²).
        # Off-diagonal elements are the covariance (ρ * σ_v * σ_f),
        # which in our case is ρ * σ² since both trait SDs are the same.
        Σ = [trait_sd^2               trait_correlation*trait_sd^2;
             trait_correlation*trait_sd^2 trait_sd^2]

        # 1c. Create the multivariate normal distribution object.
        dist = MvNormal(μ, Σ)

        # 1d. Sample all traits at once. This returns a 2xS matrix.
        traits_matrix = rand(dist, S)
        f_traits = traits_matrix[1, :] # Foraging traits (row 1)
        v_traits = traits_matrix[2, :] # Vulnerability traits (row 2)
        
        # 1e. Enforce the ecological rule that producers do not forage.
        producer_indices = findall(mc -> mc == :producer, metabolic_classes)
        f_traits[producer_indices] .= 0.0
    else
        # CASE B: Traits WERE provided. Use them directly.
        # This branch is used by `create_invaded_community`.
        v_traits = vulnerability_traits
        f_traits = foraging_traits
        # We assume the calling function (e.g., create_invaded_community)
        # has already correctly set producer foraging traits to 0.0.
    end
    # --- End of MODIFICATION 2 ---
    
    # --- BLOCK 2: Derive ltm Coefficients from Ecologically Intuitive Parameters ---
    # This block translates user-friendly inputs (like x_opt, sigma_x) into the
    # mathematical coefficients (α, β, γ) required for the model's
    # quadratic equation: log-odds = α + βX + γX² (where X is log10(mass_ratio)).
    log_x_opt = log10(x_opt)
    _gamma = -1.0 / (2.0 * sigma_x^2)
    _beta = -2.0 * _gamma * log_x_opt
    logit_P_max_target = log(P_max_target / (1.0 - P_max_target))
    _alpha = logit_P_max_target - (_beta * log_x_opt) - (_gamma * log_x_opt^2)

    # --- BLOCK 3: Main Generation Loop ---
    # This loop will run multiple times *only if* `ensure_viability` and
    # `remove_flawed` are both true, and the first attempt fails.
    max_attempts_internal = 100
    for attempt in 1:max_attempts_internal
        # --- 3a. Calculate Interaction Probabilities ---
        # Initialise the matrix to store probabilities (0.0 to 1.0).
        prob_matrix = zeros(Float64, S, S)
        
        # Iterate over every possible predator-prey pair.
        for pred_idx in species_indices
            # Producers cannot be predators, so skip them.
            if metabolic_classes[pred_idx] == :producer || bodymasses[pred_idx] <= 0; continue; end
            for prey_idx in species_indices
                # Prevent self-predation (cannibalism).
                if pred_idx == prey_idx || bodymasses[prey_idx] <= 0; continue; end
                
                # Calculate the log10 of the mass ratio (prey mass / predator mass).
                log_mass_ratio = log10(bodymasses[prey_idx] / bodymasses[pred_idx])
                if !isfinite(log_mass_ratio); continue; end # Safety check
                
                # Calculate the two components of the log-odds.
                # 1. The body-size component (from the quadratic formula).
                body_size_log_odds = _alpha + (_beta * log_mass_ratio) + (_gamma * log_mass_ratio^2)
                # 2. The latent-trait component (from trait matching).
                latent_trait_log_odds = δ * v_traits[prey_idx] * f_traits[pred_idx]
                
                # Combine them to get the total log-odds of an interaction.
                total_log_odds = body_size_log_odds + latent_trait_log_odds
                
                # Convert from log-odds back to a probability (0-1) using the logistic function.
                prob_matrix[pred_idx, prey_idx] = 1.0 / (1.0 + exp(-total_log_odds))
            end
        end

        # --- 3b. Build the Binary Adjacency Matrix (0s and 1s) ---
        binary_matrix = zeros(Int, S, S)
        if stochastic
            # STOCHASTIC MODE: Flip a weighted coin for every possible link.
            # A link (1) is formed if a random number [0,1] is less than the probability.
            for i in 1:S, j in 1:S; if rand() < prob_matrix[i, j]; binary_matrix[i, j] = 1; end; end
        else
            # DETERMINISTIC MODE: Create a web with a fixed number of links.
            # Calculate the total expected number of links by summing all probabilities.
            expected_links = round(Int, sum(filter(isfinite, prob_matrix)))
            if expected_links > 0 && isfinite(expected_links)
                # Find the indices of the `expected_links` *highest* probability pairs.
                num_to_take = min(expected_links, length(prob_matrix))
                top_indices = partialsortperm(vec(prob_matrix), 1:num_to_take, rev=true)
                # Set only those top-probability positions to 1.
                binary_matrix[top_indices] .= 1
            end
        end

        # --- 3c. Handle Viability Checks based on user flags ---
        if !ensure_viability
            # Viability checks are OFF. Return the first web immediately, regardless
            # of its structure.
            return (status=:unchecked, binary_matrix=binary_matrix, probability_matrix=prob_matrix, v_traits=v_traits, f_traits=f_traits)
        else
            # Viability checks are ON. Create a graph and test it.
            g = SimpleDiGraph(binary_matrix)
            if is_viable(g)
                # SUCCESS: The web is viable. Return it.
                return (status=:viable, binary_matrix=binary_matrix, probability_matrix=prob_matrix, v_traits=v_traits, f_traits=f_traits)
            else 
                # FAILURE: The web is not viable.
                if !remove_flawed
                    # We are told NOT to remove flawed webs. Return this
                    # flawed web immediately for analysis.
                    return (status=:structurally_flawed, binary_matrix=binary_matrix, probability_matrix=prob_matrix, v_traits=v_traits, f_traits=f_traits)
                end
                # If we ARE supposed to remove flawed webs (`remove_flawed = true`),
                # this `else` block finishes, and the `for` loop continues to the
                # next attempt to try and generate a new one.
            end
        end
    end # End of attempt loop

    # --- BLOCK 4: Handle Generation Failure ---
    # This code is only reached if `ensure_viability` and `remove_flawed`
    # were both `true`, and the loop completed `max_attempts_internal` times
    # without finding a single viable web.
    @warn "ltm_matrix failed to generate a viable food web after $(max_attempts_internal) attempts."
    return (status=:generation_failed, binary_matrix=zeros(Int,S,S), probability_matrix=zeros(Float64,S,S), v_traits=v_traits, f_traits=f_traits)
end

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
    # CHECK 1: The graph cannot be empty.
    if S == 0; return false; end

    # CHECK 2: There must be at least one producer.
    # A producer is defined by its structure: it has no outgoing links (it eats nothing).
    producers = findall(i -> outdegree(g, i) == 0, 1:S)
    if isempty(producers); return false; end

    # CHECK 3: All producers must be consumed by something.
    # This prevents producers from being dead-ends in the energy flow.
    for p in producers
        if indegree(g, p) == 0; return false; end
    end

    # CHECK 4: No species can be completely disconnected from the web.
    for i in 1:S
        if degree(g, i) == 0; return false; end
    end

    # If all checks pass, the web is structurally viable.
    return true
end

# --- 4. Standardised Wrapper Function (UPDATED) ---
"""
    generate_ltm_model(S, bodymasses, metabolic_classes)

Wrapper for the ltm function to be called from `main.jl`.
It runs the ltm generator with viability checks and then filters
the output based on the global BASAL_RANGE.
"""
function generate_ltm_model(S, bodymasses, metabolic_classes)

    # Call the ltm function, ensuring viability.
    result = ltm(1:S, bodymasses, metabolic_classes,    
                 ensure_viability=true, remove_flawed=true)
    
    # Handle internal generation failure
    if result.status != :viable
        # @warn "LTM model failed internal viability check. Skipping."
        return nothing
    end

    # Calculate emergent % basal species
    percent_basal = calculate_emergent_producers(result.binary_matrix)

    # New Requirement: Filter by emergent basal range
    if !is_in_basal_range(percent_basal)
        # @warn "LTM web failed basal range check ($percent_basal). Skipping."
        return nothing
    end

    return (adj=result.binary_matrix, percent_basal=percent_basal)
end

