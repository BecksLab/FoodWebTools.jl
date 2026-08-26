"""
    lmatrix(
        spp_list::Vector{<:Any},
        bodymass::Vector{Float64},
        is_producer::Vector{Bool};
        Ropt::Float64 = 100.0,
        γ::Float64 = 2.0,
        threshold::Float64 = 0.01,
        probabilistic::Bool = false
    )

A port of the ATNr L matrix function from Gauzens et al. (2023) based on the 
original descriptions from Schneider et al. (2016). Interactions are determined 
by allometric rules and a Ricker function defined by `Ropt` and `γ`.

Requires a species list, a vector of absolute body masses, and a boolean vector 
indicating if each species is a producer.

# Matrix Convention
- **Rows (i):** Consumers / Predators (who is doing the eating)
- **Columns (j):** Resources / Prey (who is being eaten)
- `matrix[i, j]` defines the interaction where predator `i` eats prey `j`.
- *Note:* This is the **transpose** of the default orientation in the R `ATNr` package.

# Arguments
- `probabilistic::Bool`: If `true`, returns a `Matrix{Float64}` containing raw interaction 
  probabilities. If `false` (default), returns a binary `Matrix{Int}` thresholded by `threshold`.
"""
function lmatrix(
    spp_list::Vector{<:Any},
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64 = 100.0,
    γ::Float64 = 2.0,
    threshold::Float64 = 0.01,
    probabilistic::Bool = false,
)
    # --- Initialization ---
    S = length(spp_list)
    
    # Pre-allocate based on requested output type
    if probabilistic
        link_matrix = zeros(Float64, (S, S))
    else
        link_matrix = zeros(Int, (S, S))
    end

    # --- Generate Links based on Allometric Rule ---
    for i = 1:S # Rows = Potential Predators
        # Skip if the species is a producer (cannot be a predator)
        if !is_producer[i]
            for j = 1:S # Columns = Potential Prey
                # Calculate the scaled mass ratio
                l = bodymass[i] / (bodymass[j] * Ropt)
                
                # Apply the Ricker-like function to determine link potential
                L = (l * exp(1 - l))^γ
                
                if probabilistic
                    link_matrix[i, j] = L
                else
                    if L > threshold
                        link_matrix[i, j] = 1
                    end
                end
            end 
        end 
    end 

    return link_matrix
end