"""
    adbm_parameters(bodymass::Vector{Float64}, is_producer::Vector{Bool}; kwargs...)

Initialize and validate the parameter dictionary required for the Allometric Diet 
Breadth Model (ADBM). 

# Arguments
- `bodymass::Vector{Float64}`: Absolute body masses of all species in the system.
- `is_producer::Vector{Bool}`: Boolean flags indicating whether each species is a primary producer.

# Keywords
- `e::Float64 = 1.0`: Energy content scaling coefficient.
- `a_adbm::Float64 = 0.0189`: Attack rate scaling coefficient.
- `ai::Float64 = -0.491`: Attack rate scaling exponent for prey mass (M_i).
- `aj::Float64 = -0.465`: Attack rate scaling exponent for predator mass (M_j).
- `b::Float64 = 0.401`: Critical mass ratio limit for the ratio-dependent handling time method.
- `h_adbm::Float64 = 1.0`: Handling time scaling coefficient.
- `hi::Float64 = 1.0`: Handling time scaling exponent for prey mass.
- `hj::Float64 = 1.0`: Handling time scaling exponent for predator mass.
- `n::Float64 = 1.0`: Numerical abundance scaling coefficient.
- `ni::Float64 = -0.75`: Numerical abundance scaling exponent (M^{-0.75}).
- `Hmethod::Symbol = :ratio`: Handling time formulation calculation method (`:ratio` or `:power`).
- `Nmethod::Symbol = :biomass`: Method for tracking prey abundance (`:original` allometric scaling or empirical `:biomass`).

# Returns
- `Dict{Symbol, Any}`: Clean parameter mapping ready for the ADBM engine.
"""
function adbm_parameters(
    bodymass::Vector{Float64},
    is_producer::Vector{Bool};
    e::Float64 = 1.0,      
    a_adbm::Float64 = 0.0189,
    ai::Float64 = -0.491,  
    aj::Float64 = -0.465,  
    b::Float64 = 0.401,    
    h_adbm::Float64 = 1.0,  
    hi::Float64 = 1.0,    
    hj::Float64 = 1.0,    
    n::Float64 = 1.0,      
    ni::Float64 = -0.75,   
    Hmethod::Symbol = :ratio, 
    Nmethod::Symbol = :biomass
)
    parameters = Dict{Symbol,Any}(
        :e => e, :a_adbm => a_adbm, :ai => ai, :aj => aj, :b => b,
        :h_adbm => h_adbm, :hi => hi, :hj => hj, :n => n, :ni => ni,
    )

    if Hmethod ∈ [:ratio, :power]
        parameters[:Hmethod] = Hmethod
    else
        error("Invalid value for Hmethod -- must be :ratio or :power")
    end
    
    if Nmethod ∈ [:original, :biomass]
        parameters[:Nmethod] = Nmethod
    else
        error("Invalid value for Nmethod -- must be :original or :biomass")
    end

    S = length(bodymass)
    parameters[:costMat] = ones(Float64, (S, S))
    parameters[:is_producer] = is_producer
    parameters[:bodymass] = bodymass

    return parameters
end

"""
    _get_adbm_terms(S::Int64, parameters::Dict{Symbol, Any}, biomass::Vector{Float64})

Internal helper to calculate intermediate physiological matrices: Energy content (`E`), 
Encounter rate matrix (`λ`), and Handling time matrix (`H`).

# Returns
- `Dict{Symbol, Any}`: Contains calculated arrays mapping dimensions across all S x S links.
"""
function _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    E = parameters[:e] .* parameters[:bodymass]

    if parameters[:Nmethod] == :original
        N = parameters[:n] .* (parameters[:bodymass] .^ parameters[:ni])
    elseif parameters[:Nmethod] == :biomass
        N = biomass
    end 

    # Calculate base encounter rates via outer product matrix scaling
    A_adbm = parameters[:a_adbm] *
             (parameters[:bodymass] .^ parameters[:aj]) * 
             (parameters[:bodymass] .^ parameters[:ai])'  

    for i = 1:S 
        A_adbm[:, i] = A_adbm[:, i] .* N[i] 
    end
    λ = A_adbm 

    if parameters[:Hmethod] == :ratio
        H = zeros(Float64, (S, S))
        ratios = (parameters[:bodymass] ./ parameters[:bodymass]')'
        for i = 1:S, j = 1:S 
            if ratios[j, i] < parameters[:b] 
                H[j, i] = parameters[:h_adbm] / (parameters[:b] - ratios[j, i])
            else
                H[j, i] = Inf
            end
        end
    elseif parameters[:Hmethod] == :power
        H = parameters[:h_adbm] *
            (parameters[:bodymass] .^ parameters[:hj]) * 
            (parameters[:bodymass] .^ parameters[:hi])'  
    end 

    return Dict{Symbol,Any}(:E => E, :λ => λ, :H => H)
end 

"""
    _get_feeding_links(E, λ, H, biomass, j)

Calculate optimal diet choices for a single targeted predator index `j` using 
the steps of optimal foraging theory. 

Filters out links tracking extinct prey items showing zero current biomass.
"""
function _get_feeding_links(
    E::Vector{Float64}, 
    λ::Matrix{Float64}, 
    H::Matrix{Float64}, 
    biomass::Vector{Float64}, 
    j::Int,             
)
    profit = E ./ H[j, :]
    profit[biomass .== 0.0] .= -1.0 

    profs = sortperm(profit, rev = true)

    λSort = λ[j, profs] 
    HSort = H[j, profs] 
    ESort = E[profs]    

    λH = cumsum(λSort .* HSort)
    Eλ = cumsum(ESort .* λSort)

    λH[isnan.(λH)] .= Inf
    Eλ[isnan.(Eλ)] .= Inf

    cumulativeProfit = Eλ ./ (1 .+ λH)

    if all(0 .== cumulativeProfit) || all(isnan.(cumulativeProfit))
        return Int[]
    else
        max_profit_indices = findall(cumulativeProfit .== maximum(filter(!isnan, cumulativeProfit)))
        optimal_breadth = maximum(max_profit_indices)
        feeding = profs[1:optimal_breadth]
        
        # Explicitly remove any prey that have zero biomass 
        # to prevent creating structural links to extinct nodes.
        return filter(p -> biomass[p] > 0.0, feeding)
    end
end 

"""
    adbm(spp_list::Vector{Any}, parameters::Dict{Symbol, Any}, biomass::Vector{Float64})

Generate a structural food web adjacency matrix matching the Allometric Diet Breadth Model 
framework described by Petchey et al. (2008).

# Matrix Structure (CRITICAL)
- **Rows (i):** Consumers / Predators (who is doing the eating).
- **Columns (j):** Resources / Prey (who is being eaten).
- `adbmMAT[i, j] = true` indicates predator `i` consumes prey `j`.

# Returns
- `Matrix{Bool}`: S x S binary topology web network.
"""
function adbm(spp_list::Vector{Any}, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    S = length(spp_list)
    adbmMAT = zeros(Bool, (S, S))

    adbmTerms = _get_adbm_terms(S, parameters, biomass)
    E = adbmTerms[:E]
    λ = adbmTerms[:λ]
    H = adbmTerms[:H]

    for j = 1:S
        if !parameters[:is_producer][j] && biomass[j] > 0.0
            feeding = _get_feeding_links(E, λ, H, biomass, j)
            if !isempty(feeding)
                adbmMAT[j, feeding] .= true 
            end
        end
    end 

    return adbmMAT
end
