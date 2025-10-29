"""
adbm.jl
-----------------
Generates a food web using the Allometric Diet Breadth Model (ADBM)
from Petchey et al. (2008), using functions ported from
BioEnergeticFoodWebs.jl.

Contains the core ADBM helper functions (`adbm_parameters`, `_get_adbm_terms`,
`_get_feeding_links`), the main `adbm` generation function (provided by user),
and a standardised wrapper function `generate_adbm_model` called by `main.jl`.

The wrapper ensures the correct biomass vector (`M^-0.75`) is used via
`Nmethod = :biomass` and applies external checks for basal percentage and
connectance.
"""

# --- 1. Dependencies ---
using Statistics # Needed for cumsum

# --- 2. ADBM Helper Functions (User Provided Code) ---

"""
  adbm_parameters()

  Returns the parameters needed for the adbm model. Defaults to the values
  specified in BioEnergeticFoodWebs.jl.
"""
function adbm_parameters(
    bodymass::Vector{Float64}, # Expects absolute body masses
    is_producer::Vector{Bool};
    # Default parameters based on BioEnergeticFoodWebs.jl / common usage
    e::Float64 = 1.0,      # Energy content scaling (often assumed constant)
    a_adbm::Float64 = 0.0189,# Attack rate coefficient
    ai::Float64 = -0.491,  # Attack rate scaling exponent for prey mass
    aj::Float64 = -0.465,  # Attack rate scaling exponent for predator mass
    b::Float64 = 0.401,    # Critical mass ratio for handling time (ratio method)
    h_adbm::Float64 = 1.0,  # Handling time coefficient
    hi::Float64 = 1.0,    # Handling time scaling exponent for prey mass (power method)
    hj::Float64 = 1.0,    # Handling time scaling exponent for predator mass (power method)
    n::Float64 = 1.0,      # Abundance coefficient (original method)
    ni::Float64 = -0.75,   # Abundance scaling exponent (original method, M^-0.75)
    Hmethod::Symbol = :ratio, # Method for calculating handling time
    Nmethod::Symbol = :original # Method for calculating prey availability (N)
)
    # --- Initialize Parameter Dictionary ---
    parameters = Dict{Symbol,Any}(
        :e => e, :a_adbm => a_adbm, :ai => ai, :aj => aj, :b => b,
        :h_adbm => h_adbm, :hi => hi, :hj => hj, :n => n, :ni => ni,
    )

    # --- Validate Method Keywords ---
    # Check if the chosen handling time method is valid.
    if Hmethod ∈ [:ratio, :power]
        parameters[:Hmethod] = Hmethod
    else
        error("Invalid value for Hmethod -- must be :ratio or :power")
    end
    # Check if the chosen prey availability method is valid.
    if Nmethod ∈ [:original, :biomass]
        parameters[:Nmethod] = Nmethod
    else
        error("Invalid value for Nmethod -- must be :original or :biomass")
    end

    # --- Add Other Necessary Parameters ---
    S = length(bodymass)
    # Add a cost matrix (usually assumed to be all 1s if not specified).
    parameters[:costMat] = ones(Float64, (S, S))
    # Store the producer boolean vector.
    parameters[:is_producer] = is_producer
    # Store the absolute body mass vector.
    parameters[:bodymass] = bodymass

    return parameters
end # end adbm_parameters

"""
  _get_adbm_terms(...)
Calculates intermediate terms needed for the ADBM optimal foraging calculation,
namely energy content (E), encounter rate matrix (λ), and handling time
matrix (H), based on the provided parameters and biomass/abundance method.
... (full original docstring) ...
"""
function _get_adbm_terms(S::Int64, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})
    # --- Calculate Energy Content (E) ---
    # Assumes energy content scales linearly with body mass (e=1 default).
    E = parameters[:e] .* parameters[:bodymass]

    # --- Calculate Prey Availability (N) ---
    # Based on the chosen Nmethod.
    if parameters[:Nmethod] == :original
        # Original method: N scales allometrically with prey body mass.
        N = parameters[:n] .* (parameters[:bodymass] .^ parameters[:ni])
    elseif parameters[:Nmethod] == :biomass
        # Biomass method: N is directly the provided biomass vector.
        N = biomass
    end # Nmethod check ignored as wrapper forces :biomass

    # --- Calculate Encounter Rate Base (A_adbm) ---
    # Calculates the base encounter rate without prey availability scaling.
    # Uses matrix operations for efficiency: (pred_mass^aj) * (prey_mass^ai)'
    A_adbm =
        parameters[:a_adbm] *
        (parameters[:bodymass] .^ parameters[:aj]) * # Predator scaling (Column vector)
        (parameters[:bodymass] .^ parameters[:ai])'  # Prey scaling (Row vector)

    # --- Calculate Full Encounter Rate Matrix (λ) ---
    # Scale the base rate by prey availability (N).
    for i = 1:S # Iterate through prey (columns)
        A_adbm[:, i] = A_adbm[:, i] .* N[i] # Multiply predator column by prey availability
    end
    λ = A_adbm # λ[pred, prey]

    # --- Calculate Handling Time Matrix (H) ---
    # Based on the chosen Hmethod.
    if parameters[:Hmethod] == :ratio
        # --- Ratio Method ---
        H = zeros(Float64, (S, S))
        # Calculate matrix of PreyMass / PredMass ratios.
        # Transpose ensures Preds in Rows, Prey in Cols to match H structure.
        ratios = (parameters[:bodymass] ./ parameters[:bodymass]')'
        for i = 1:S, j = 1:S # i = prey, j = predator
            # Handling time increases dramatically as prey size approaches critical ratio 'b'.
            if ratios[j, i] < parameters[:b] # If Prey/Pred ratio < b
                H[j, i] = parameters[:h_adbm] / (parameters[:b] - ratios[j, i])
            else
                # Handling time is infinite if prey is too large relative to predator.
                H[j, i] = Inf
            end
        end
    elseif parameters[:Hmethod] == :power
        # --- Power Method ---
        # Handling time scales allometrically with predator and prey mass.
        H =
            parameters[:h_adbm] *
            (parameters[:bodymass] .^ parameters[:hj]) * # Predator scaling
            (parameters[:bodymass] .^ parameters[:hi])'  # Prey scaling
    end # end Hmethod check

    # --- Store and Return Terms ---
    adbmTerms = Dict{Symbol,Any}(:E => E, :λ => λ, :H => H)
    return adbmTerms
end # end _get_adbm_terms

"""
  _get_feeding_links(...)
Determines the optimal diet for a single predator `j` based on the
profitability of available prey, encounter rates, and handling times,
following the marginal value theorem logic.
... (full original docstring) ...
"""
function _get_feeding_links(
    E::Vector{Float64}, # Energy content vector for all species
    λ::Matrix{Float64}, # Encounter rate matrix [pred, prey]
    H::Matrix{Float64}, # Handling time matrix [pred, prey]
    biomass::Vector{Float64}, # Biomass vector (used to exclude extinct prey)
    j::Int,             # Index of the predator species
)
    # --- Calculate Profitability ---
    # Profitability = Energy Gain / Handling Time for predator j eating each prey.
    profit = E ./ H[j, :] # Takes the j-th row of H
    # Exclude prey with zero biomass by setting their profitability to -1.
    profit[vec(biomass .== 0.0)] .= -1.0

    # --- Rank Prey by Profitability ---
    # Get the indices that would sort the `profit` vector in descending order.
    profs = sortperm(profit, rev = true)

    # --- Sort E, λ, H according to Profitability Rank ---
    # Reorder encounter rates, handling times, and energy content for predator j
    # based on the profitability ranking.
    λSort = λ[j, profs] # Encounter rates sorted by profitability
    HSort = H[j, profs] # Handling times sorted by profitability
    ESort = E[profs]    # Energy contents sorted by profitability

    # --- Calculate Cumulative Terms ---
    # Calculate cumulative sum of (Encounter Rate * Handling Time).
    λH = cumsum(λSort .* HSort)
    # Calculate cumulative sum of (Energy Content * Encounter Rate).
    Eλ = cumsum(ESort .* λSort)

    # Handle potential NaNs or Infs resulting from division by zero or Inf handling times.
    λH[isnan.(λH)] .= Inf
    Eλ[isnan.(Eλ)] .= Inf

    # --- Calculate Cumulative Profit Rate ---
    # This represents the overall energy intake rate if the diet includes
    # prey up to the current rank (based on marginal value theorem logic).
    # Rate = Cumulative Energy Gain / (1 + Cumulative Handling Time * Encounter Rate)
    cumulativeProfit = Eλ ./ (1 .+ λH)

    # --- Determine Optimal Diet Breadth ---
    # Find the rank(s) that yield the maximum cumulative profit rate.
    if all(0 .== cumulativeProfit) || all(isnan.(cumulativeProfit))
        # If all profit rates are zero or NaN, the optimal diet is empty.
        feeding = []
    else
        # Find all indices where the maximum profit rate occurs.
        max_profit_indices = findall(cumulativeProfit .== maximum(filter(!isnan, cumulativeProfit)))
        # The optimal diet includes all prey up to the *last* index that yields
        # the maximum profit rate.
        optimal_breadth = maximum(max_profit_indices)
        # Select the original indices (from `profs`) corresponding to this breadth.
        feeding = profs[1:optimal_breadth]
    end

    # --- Return Optimal Prey Indices ---
    return feeding # Vector of indices of prey included in the optimal diet
end # end _get_feeding_links

# --- 3. ADBM Main Generator (User Provided Code) ---
"""
adbmmodel(spp_list::Vector{Any}, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})

  This function returns the food web based on the ADBM model of Petchey et al. 2008.
  The function takes the paramteres created by adbm_parameters and uses 
  _get_adbm_terms and _get_feeding_links to determine the web structure.

  Note this (and all internal) functions has been ported from the 
  BioEnergeticFoodWebs.jl source code and has been (minimally) modified for the 
  purpose of this project.

  #### References

  Petchey, Owen L., Andrew P. Beckerman, Jens O. Riede, and Philip H. Warren.
  2008. “Size, Foraging, and Food Web Structure.” Proceedings of the National
  Academy of Sciences 105 (11): 4191–96. https://doi.org/10.1073/pnas.0710672105.

"""
function adbm(spp_list::Vector{Any}, parameters::Dict{Symbol,Any}, biomass::Vector{Float64})

    # --- Initialization ---
    S = length(spp_list)
    # Initialize Boolean adjacency matrix.
    adbmMAT = zeros(Bool, (S, S))

    # --- Pre-calculate ADBM Terms ---
    # Calculate E, λ, H matrices once.
    adbmTerms = _get_adbm_terms(S, parameters, biomass)
    E = adbmTerms[:E]
    λ = adbmTerms[:λ]
    H = adbmTerms[:H]

    # --- Determine Feeding Links for Each Predator ---
    # Iterate through each species index `j`.
    for j = 1:S
        # Skip if the species is a producer.
        if !parameters[:is_producer][j]
            # Also skip if the predator has zero biomass (cannot forage).
            if biomass[j] > 0.0
                # Calculate the optimal diet (indices of prey) for predator j.
                feeding = _get_feeding_links(E, λ, H, biomass, j)
                # If the diet is not empty, set corresponding entries in the matrix row.
                if !isempty(feeding)
                    adbmMAT[j, feeding] .= 1 # Predator j eats Prey in `feeding`
                end
            end
        end
    end # end predator loop

    # --- Return Boolean Matrix ---
    return adbmMAT
end # end adbm

