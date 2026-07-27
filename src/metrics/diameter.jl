"""
    diameter(A::Matrix{Bool}) -> Int

Calculate the diameter of a directed food web.

The input is interpreted as a directed adjacency matrix, where

    A[i, j] = true

indicates an edge from node `i` to node `j`.

The diameter is defined as the length of the longest finite shortest path
between any pair of vertices. Unreachable vertex pairs are ignored when
computing the maximum.

Shortest paths are computed by performing a breadth-first search (BFS)
from every vertex.

# Arguments

- `A::Matrix{Bool}`: directed adjacency matrix.

# Returns

The diameter of the graph as an integer number of edges.

Returns `0` if the graph contains no edges or no reachable paths between
distinct vertices.

# Examples

```julia
A = Bool[
    0 1 0 0
    0 0 1 0
    0 0 0 1
    0 0 0 0
]

diameter(A) == 3
```
"""
function diameter(A::AbstractMatrix{Bool})
    n = size(A, 1)
    dist = fill(-1, n)
    queue = Vector{Int}(undef, n)

    diam = 0

    for src in 1:n
        fill!(dist, -1)
        dist[src] = 0

        front = back = 1
        queue[1] = src

        while front <= back
            u = queue[front]
            front += 1

            @inbounds for v in 1:n
                if A[u, v] && dist[v] == -1
                    dist[v] = dist[u] + 1
                    back += 1
                    queue[back] = v
                    diam = max(diam, dist[v])
                end
            end
        end
    end

    return diam
end