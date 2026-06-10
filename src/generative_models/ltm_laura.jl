# ========================================================================= #
#    LATENT TRAIT MODEL (LAURA'S CONFIGURATION)                             #
# ========================================================================= #

using Distributions, Random

"""
    is_mass_ratio_valid(pred_mass::Float64, prey_mass::Float64, envelope::Tuple{Float64, Float64}) -> Bool

Acts as a 'Topological Gatekeeper' to prevent the generation of physically impossible 
trophic links. It determines if a predator-prey mass ratio (R = predator_mass / prey_mass) 
is biologically feasible within a user-defined physical envelope.

# Arguments
- `pred_mass::Float64`: The body mass of the attacking predator (g).
- `prey_mass::Float64`: The body mass of the target prey organism (g).
- `envelope::Tuple{Float64, Float64}`: A tuple representing the [Lower, Upper] bounds 
  of the allowed mass ratio. 
"""
function is_mass_ratio_valid(pred_mass::Float64, prey_mass::Float64, envelope::Tuple{Float64, Float64})
    # --- Biological Sanity Check ---
    if pred_mass <= 0.0 || prey_mass <= 0.0
        return false
    end

    # --- Calculate Trophic Mass Ratio (R) ---
    ratio = pred_mass / prey_mass
    
    # --- Enforce the Physical Envelope ---
    if ratio >= envelope[1] && ratio <= envelope[2]
        return true
    else
        return false
    end
end

"""
    ltm_laura(...) -> NamedTuple

Generates a food web interaction matrix using a Latent Trait Model (LTM) subject 
to strict allometric bounding. Designed to slot directly into established testing 
pipelines by mirroring the return structure of the original `ltm` function.
"""
function ltm_laura(
    species_indices::AbstractVector{Int},
    bodymasses::AbstractVector{Float64},
    metabolic_classes::AbstractVector{Symbol}; 
    
    # LTM Parameters
    P_max_target::Float64 = 0.85,
    x_opt::Float64 = 100.0, 
    sigma_x::Float64 = 1.5,
    δ::Float64 = 8.0,
    
    # Physical Constraints
    mass_ratio_envelope::Tuple{Float64, Float64} = (1e-1, 1e4),
    
    # Trait Generation
    trait_sd::Float64 = 1.0,
    trait_correlation::Float64 = 1.0, 
    vulnerability_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    foraging_traits::Union{AbstractVector{Float64}, Nothing} = nothing,
    
    # Generation Settings
    stochastic::Bool = true,       
    ensure_viability::Bool = false,  
    remove_flawed::Bool = false
)
    # --- Initialization ---
    S = length(species_indices)
    local f_traits, v_traits 

    # --- BLOCK 1: Latent Trait Generation (Dynamic Routing) ---
    if isnothing(vulnerability_traits) || isnothing(foraging_traits)
        
        # ROUTE A: Absolute 100% Correlation (Mathematically Safe Path)
        if trait_correlation == 1.0
            dist = Normal(0.0, trait_sd)
            base_traits = rand(dist, S) 
            v_traits = copy(base_traits) 
            f_traits = copy(base_traits) 
            
        # ROUTE B: Partial Correlation (Multivariate Normal Path)
        else
            μ = [0.0, 0.0]
            Σ = [trait_sd^2               trait_correlation*trait_sd^2;
                 trait_correlation*trait_sd^2 trait_sd^2]
            dist = MvNormal(μ, Σ)
            traits_matrix = rand(dist, S) 
            f_traits = traits_matrix[1, :] 
            v_traits = traits_matrix[2, :] 
        end
        
        # Override: Primary producers do not forage.
        producer_indices = findall(mc -> mc == :producer, metabolic_classes)
        f_traits[producer_indices] .= 0.0
        
    else
        # Use Provided Traits
        v_traits = vulnerability_traits
        f_traits = foraging_traits
    end 

    # --- BLOCK 2: Derive LTM Coefficients ---
    log_x_opt = log10(x_opt)
    _gamma = -1.0 / (2.0 * sigma_x^2)
    _beta = -2.0 * _gamma * log_x_opt
    logit_P_max_target = log(P_max_target / (1.0 - P_max_target))
    _alpha = logit_P_max_target - (_beta * log_x_opt) - (_gamma * log_x_opt^2)

    # --- BLOCK 3: Generate Interaction Probabilities ---
    prob_matrix = zeros(Float64, S, S)
    
    for pred_idx in species_indices
        if metabolic_classes[pred_idx] == :producer || bodymasses[pred_idx] <= 0; continue; end
        
        for prey_idx in species_indices
            if pred_idx == prey_idx || bodymasses[prey_idx] <= 0; continue; end
            
            # ✨ TOPOLOGICAL GATEKEEPER INJECTION ✨
            if !is_mass_ratio_valid(bodymasses[pred_idx], bodymasses[prey_idx], mass_ratio_envelope)
                prob_matrix[pred_idx, prey_idx] = 0.0
                continue 
            end
            
            # Calculate log-odds
            log_mass_ratio = log10(bodymasses[pred_idx] / bodymasses[prey_idx])
            if !isfinite(log_mass_ratio); continue; end 
            
            body_size_log_odds = _alpha + (_beta * log_mass_ratio) + (_gamma * log_mass_ratio^2)
            latent_trait_log_odds = δ * v_traits[prey_idx] * f_traits[pred_idx]
            total_log_odds = body_size_log_odds + latent_trait_log_odds
            
            # Convert to probability
            prob_matrix[pred_idx, prey_idx] = 1.0 / (1.0 + exp(-total_log_odds))
        end 
    end 

    # --- BLOCK 4: Build Binary Adjacency Matrix ---
    binary_matrix = zeros(Int, S, S)
    
    if stochastic 
        # Stochastic Realisation
        for i in 1:S, j in 1:S
            if rand() < prob_matrix[i, j]; binary_matrix[i, j] = 1; end
        end
    else
        # Deterministic Realisation
        expected_links = round(Int, sum(filter(isfinite, prob_matrix)))
        if expected_links > 0 && isfinite(expected_links)
            finite_indices = findall(isfinite, vec(prob_matrix))
            num_finite = length(finite_indices)
            num_to_take = min(expected_links, num_finite)
            
            if num_to_take > 0
                finite_probs = vec(prob_matrix)[finite_indices]
                p = sortperm(finite_probs, rev=true)
                top_indices_relative = p[1:num_to_take]
                top_indices = finite_indices[top_indices_relative]
                binary_matrix[top_indices] .= 1
            end 
        end 
    end 

    # --- BLOCK 5: Return Result ---
    # Bypassing strict graph evaluation to match original ltm integration structure
    return (status=:unchecked, binary_matrix=binary_matrix, probability_matrix=prob_matrix, v_traits=v_traits, f_traits=f_traits)

end