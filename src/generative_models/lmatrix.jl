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
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    Ropt::Float64 = 100.0,
    γ::Float64 = 2.0,
    threshold::Float64 = 0.01,
)

    S = length(spp_list)

    prob_matrix = zeros(Bool, (S, S))

    for i = 1:S
        # only assess if consumer is not a producer
        if !is_producer[i]
            for j = 1:S
                l = bodymass[i] / (bodymass[j] * Ropt)
                L = (l * exp(1 - l))^γ
                if L > threshold
                    prob_matrix[i, j] = 1
                end
            end
        end
    end

    return prob_matrix
end
