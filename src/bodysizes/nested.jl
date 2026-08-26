# for dimulating bodysizes following the appraoch of Stouffer and Bascompte 2010
# doi: /10.1111/j.1461-0248.2009.01407.x

using Distributions

"""
    initial_bodymasses(classes, size_bounds, dist)

Draw one body mass per species by sampling from `dist` truncated to the
bounds associated with each class.

# Arguments
- `classes`: Vector of size-class labels (e.g. `["small","large","small"]`).
- `size_bounds`: Dictionary mapping class labels to `(min, max)` bounds.
- `dist`: Distribution from which body masses are drawn before truncation.

# Returns
A vector of body masses, one per species.
"""
function initial_bodymasses(classes, size_bounds, dist)
    map(classes) do cls
        lo, hi = size_bounds[cls]
        rand(truncated(dist, lo, hi))
    end
end

"""
    loglikelihood(logM, A; μ=6.1, σ=5.75)

Log-likelihood of body masses under the Stouffer & Bascompte (2010)
predator-prey mass ratio model.
"""
function loglikelihood(logM, A; μ = 6.1, σ = 5.75)
    model = Normal(μ, σ)
    ll = 0.0

    for predator in axes(A, 1), prey in axes(A, 2)
        if A[predator, prey] == 1
            ratio = logM[predator] - logM[prey]
            ll += logpdf(model, ratio)
        end
    end

    ll
end

function propose(logM, idx, bounds; proposal_sd = 0.25)
    proposal = copy(logM)

    lo, hi = log.(bounds)
    proposal[idx] = clamp(
        proposal[idx] + rand(Normal(0, proposal_sd)),
        lo,
        hi
    )

    proposal
end

function metropolis_step(logM, A, species_bounds;
                         μ = 6.1,
                         σ = 5.75,
                         proposal_sd = 0.25)

    current = loglikelihood(logM, A; μ, σ)

    idx = rand(eachindex(logM))
    proposal = propose(logM, idx, species_bounds[idx];
                       proposal_sd)

    proposed = loglikelihood(proposal, A; μ, σ)

    log(rand()) < proposed - current ? proposal : logM
end

"""
    nested_bodymasses(classes, A, size_bounds, dist; μ=6.1, σ=5.75, iterations=20_000)

Estimate nested body masses using the Metropolis-Hastings algorithm described in
Stouffer & Bascompte (2010).

The algorithm:

1. Samples an initial body mass for each species from `dist`, truncated to its
   allowed size interval.
2. Converts body masses to log space.
3. Repeatedly proposes a perturbation to one species.
4. Accepts or rejects the proposal according to the predator-prey body-mass
   ratio likelihood.

# Arguments
- `classes`: Vector assigning each species to a size class.
- `A`: Binary adjacency matrix (`A[i,j] = 1` if species `i` consumes species `j`).
- `size_bounds`: Mapping from class labels to `(lower, upper)` body mass bounds.
- `dist`: Global body-mass distribution used for initialization.

# Keyword arguments
- `μ`: Mean log predator-prey mass ratio.
- `σ`: Standard deviation of the log ratio.
- `iterations`: Number of Metropolis-Hastings iterations.

# Returns
A vector of estimated body masses.
"""
function nested_bodymasses(classes, A, size_bounds, dist;
                           μ = 6.1,
                           σ = 5.75,
                           iterations = 20_000)

    M = initial_bodymasses(classes, size_bounds, dist)
    logM = log.(M)

    species_bounds = [size_bounds[c] for c in classes]

    for _ in 1:iterations
        logM = metropolis_step(logM, A, species_bounds; μ, σ)
    end

    exp.(logM)
end