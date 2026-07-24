"""
    clustering(A::Matrix{Bool}; include_isolates=false) -> Float64

Calculate the mean local clustering coefficient of a graph.

The input is interpreted as a possibly directed adjacency matrix.
The graph is converted to an undirected graph using:

    A_undir = A ∪ A'

Self-loops are removed.

The local clustering coefficient of node `i` is:

    C_i = T_i / (k_i * (k_i - 1) / 2)

where:

- `k_i` is the number of neighbours of node `i`
- `T_i` is the number of triangles containing node `i`

Nodes with degree 0 or 1 have undefined clustering coefficients.

By default (`include_isolates=false`), these nodes are excluded from the
mean.

If `include_isolates=true`, they contribute a clustering coefficient of 0.

# Arguments

- `A::Matrix{Bool}`: adjacency matrix
- `include_isolates::Bool=false`: whether to include nodes with degree ≤ 1
  in the mean calculation

# Returns

The mean local clustering coefficient.

# Examples

```julia
A = Bool[
    0 1 1
    1 0 1
    1 1 0
]

clustering(A) == 1.0
"""
function clustering(A::Matrix{Bool}; include_isolates=false)

    N = size(A, 1)

    A_undir = A .| A'

    # remove self loops
    A_undir[diagind(A_undir)] .= false

    # degree
    k = sum(A_undir, dims=2)[:]

    # triangles involving each node
    triangle_count = diag(A_undir^3) ./ 2

    C_values = Float64[]

    for i in 1:N

        denominator = k[i] * (k[i]-1) / 2

        if denominator == 0

            if include_isolates
                push!(C_values, 0.0)
            end

            continue
        end

        push!(
            C_values,
            triangle_count[i] / denominator
        )

    end

    # If the graph has no nodes with defined clustering
    isempty(C_values) && return 0.0

    return mean(C_values)
end