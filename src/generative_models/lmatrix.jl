"""
lmatrix.jl
-----------------
Generates a food web using the Allometric Trophic Network (ATN)
L-matrix model from Gauzens et al. (2023).

Contains the core `lmatrix` generation function (provided by user) and a
standardised wrapper function `generate_lmatrix` called by `main.jl`.

Applies external checks for basal percentage and connectance in the wrapper.
"""

# --- 1. Dependencies (loaded by main.jl) ---
# (None specific to this file)

# --- 2. ATN Generator (User Provided Code) ---
"""
  lmatrix(
    spp_list::Vector{Any},
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64,
    γ::Float64,
    threshold::Float64,
)

    A port of the ATNr L matrix function from Gauzens et al. (2023) based
    on the original descriptions from Schneider et al. (2016). 
    Interactions are determined by allometric rules and a Ricker function 
    defined by `Ropt` and `γ` and returns a probabilistic 
    `SpeciesInteractionNetwork`. `Ropt`, `γ`, and `threshold` use the ATNr
    defaults.

    Requires a species list, vector of bodymass, and if the species is a
    producer or not.
    
    #### References
    
    Gauzens, B., Brose U., Delmas E., and Berti E. 2023. “ATNr: Allometric 
    Trophic Network Models in R.” Methods in Ecology and Evolution 14 (11): 
    2766–73. https://doi.org/10.1111/2041-210X.14212.

    Schneider, Florian D., Ulrich Brose, Björn C. Rall, and Christian Guill.
    2016. “Animal Diversity and Ecosystem Functioning in Dynamic Food Webs.”
    Nature Communications 7 (1): 12718. https://doi.org/10.1038/ncomms12718.

"""
function lmatrix(
    spp_list::Vector{Any},
    bodymass::Vector{Float64}, # Expects absolute body masses
    is_producer::Vector{Bool};
    Ropt::Float64 = 100.0,
    γ::Float64 = 2.0,
    threshold::Float64 = 0.01,
)
    # --- Initialization ---
    S = length(spp_list)
    # Initialize the matrix as Boolean (will be converted to Int later).
    link_matrix = zeros(Bool, (S, S))

    # --- Generate Links based on Allometric Rule ---
    # Iterate through potential predators.
    for i = 1:S
        # Skip if the species is a producer (cannot be a predator).
        if !is_producer[i]
            # Iterate through potential prey.
            for j = 1:S
                # Calculate the scaled mass ratio (predator mass / (prey mass * Ropt)).
                l = bodymass[i] / (bodymass[j] * Ropt)
                # Apply the Ricker-like function to determine link potential.
                # Note: `exp(1-l)` ensures peak at l=1 (pred_mass = prey_mass * Ropt).
                L = (l * exp(1 - l))^γ
                # Add a link if the potential exceeds the threshold.
                if L > threshold
                    link_matrix[i, j] = 1 # Predator i eats Prey j
                end
            end # end prey loop
        end # end producer check
    end # end predator loop

    # --- Return Boolean Matrix ---
    return link_matrix
end # end l_matrix

# --- 3. Standardised Wrapper Function ---
"""
    generate_lmatrix(S, bodymasses, is_producer_bool)

Wrapper for the `l_matrix` function called by `main.jl`.
Runs the ATN generator and applies external checks for basal percentage
and connectance ranges.

# Arguments
- `S::Int`: Species richness.
- `bodymasses::Vector{Float64}`: Absolute body masses.
- `is_producer_bool::Vector{Bool}`: Boolean vector indicating producers.

# Returns
- `NamedTuple`: Contains `adj`, `percent_basal`, `connectance` if checks pass.
- `nothing`: If the basal or connectance check fails.
"""
function generate_lmatrix(S, bodymasses, is_producer_bool)

    # --- Define Species List ---
    # Simple integer list from 1 to S.
    spp_list = 1:S
    
    # --- Call Core ATN Function ---
    # Explicitly pass Ropt=100 (matches default but ensures consistency).
    adj_bool = l_matrix(spp_list, bodymasses, is_producer_bool; Ropt = 100.0)
    
    # --- Convert to Integer Matrix ---
    # Convert the Boolean matrix to an Integer matrix for consistency.
    adj_int = Matrix{Int}(adj_bool)

    # --- CHECK 1: Emergent Basal % ---
    percent_basal = calculate_emergent_producers(adj_int)
    if !is_in_basal_range(percent_basal)
        return nothing # Fail if outside range.
    end

    # --- CHECK 2: Emergent Connectance ---
    connectance = calculate_connectance(adj_int)
    if !is_in_connectance_range(connectance)
        return nothing # Fail if outside range.
    end

    # --- Return Success ---
    # If both checks passed, return the results.
    return (adj=adj_int, percent_basal=percent_basal, connectance=connectance)
end # end generate_lmatrix
