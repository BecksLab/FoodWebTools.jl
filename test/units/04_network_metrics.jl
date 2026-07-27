using FoodWebTools
using Test
using LinearAlgebra
using Statistics


@testset "Clustering Tests" begin


    @testset "Complete triangle" begin

        A = Bool[
            0 1 1
            1 0 1
            1 1 0
        ]

        @test clustering(A) ≈ 1.0
        @test clustering(A; include_isolates=true) ≈ 1.0

    end


    @testset "Linear chain has zero clustering" begin

        A = Bool[
            0 1 0
            1 0 1
            0 1 0
        ]

        @test clustering(A) ≈ 0.0
        @test clustering(A; include_isolates=true) ≈ 0.0

    end


    @testset "Square lattice without diagonals" begin

        A = Bool[
            0 1 0 1
            1 0 1 0
            0 1 0 1
            1 0 1 0
        ]

        @test clustering(A) ≈ 0.0
        @test clustering(A; include_isolates=true) ≈ 0.0

    end


    @testset "Directed cycle becomes undirected triangle" begin

        A = Bool[
            0 1 0
            0 0 1
            1 0 0
        ]

        # A directed cycle becomes a triangle after A ∪ A'
        @test clustering(A) ≈ 1.0

    end


    @testset "Triangle plus isolated node" begin

        A = Bool[
            0 1 1 0
            1 0 1 0
            1 1 0 0
            0 0 0 0
        ]

        # Default behaviour: ignore isolated nodes
        @test clustering(A) ≈ 1.0

        # Including isolated nodes adds a zero contribution
        @test clustering(A; include_isolates=true) ≈ 0.75

    end


    @testset "Triangle plus node" begin

        A = Bool[
            0 1 1 0
            1 0 1 0
            1 1 0 1
            0 0 1 0
        ]

        # Degree-one node excluded, but node 3 is affected by the extra neighbour
        @test clustering(A) ≈ 7/9

        # Degree-one node included as zero
        @test clustering(A; include_isolates=true) ≈ 7/12

    end


    @testset "Self loops are ignored" begin

        A = Bool[
            1 1 1
            1 1 1
            1 1 1
        ]

        B = Bool[
            0 1 1
            1 0 1
            1 1 0
        ]

        @test clustering(A) ≈ clustering(B)

    end


    @testset "All isolated nodes" begin

        A = Bool[
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0
            0 0 0 0 0
        ]

        # No nodes have defined clustering
        @test clustering(A) ≈ 0.0

        # Including isolates gives all zeros
        @test clustering(A; include_isolates=true) ≈ 0.0

    end


    @testset "Single node graph" begin

        A = Bool[0;;]

        @test clustering(A) ≈ 0.0
        @test clustering(A; include_isolates=true) ≈ 0.0

    end


    @testset "Two node graph" begin

        A = Bool[
            0 1
            1 0
        ]

        @test clustering(A) ≈ 0.0
        @test clustering(A; include_isolates=true) ≈ 0.0

    end


    @testset "Random graphs produce valid clustering values" begin

        for n in 2:25

            A = rand(Bool, n, n)

            # remove self loops
            A[diagind(A)] .= false

            C1 = clustering(A)
            C2 = clustering(A; include_isolates=true)

            @test 0.0 ≤ C1 ≤ 1.0
            @test 0.0 ≤ C2 ≤ 1.0

        end

    end


end

@testset "Intervality Tests" begin

    @testset "swap_matrix_rows_cols!" begin
        A = Bool[
            0 1 0;
            1 0 1;
            0 1 0
        ]

        B = copy(A)
        FoodWebTools.swap_matrix_rows_cols!(B, 1, 3)

        @test B == Bool[
            0 1 0;
            1 0 1;
            0 1 0
        ]  # symmetric under swapping 1 and 3

        # Swapping a matrix twice should recover the original
        C = Bool[
            0 1 1;
            1 0 0;
            1 0 0
        ]

        original = copy(C)
        FoodWebTools.swap_matrix_rows_cols!(C, 1, 2)
        FoodWebTools.swap_matrix_rows_cols!(C, 1, 2)

        @test C == original
    end

    @testset "calculate_gaps" begin

        # No feeding links
        A = falses(4, 4)
        @test FoodWebTools.calculate_gaps(A) == 0

        # Single prey
        A = Bool[
            1 0 0;
            0 1 0;
            0 0 1
        ]
        @test FoodWebTools.calculate_gaps(A) == 0

        # Continuous diet
        A = Bool[
            1 1 1 0;
            0 1 1 0;
            0 0 1 1;
            0 0 0 0
        ]
        @test FoodWebTools.calculate_gaps(A) == 0

        # One gap
        A = Bool[
            1 0 1 0;
            0 0 0 0;
            0 0 0 0;
            0 0 0 0
        ]
        @test FoodWebTools.calculate_gaps(A) == 1

        # Two gaps
        A = Bool[
            1 0 0 1 0;
            0 0 0 0 0;
            0 0 0 0 0;
            0 0 0 0 0;
            0 0 0 0 0
        ]
        @test FoodWebTools.calculate_gaps(A) == 2

        # Multiple consumers
        A = Bool[
            1 0 1 0;  # 1 gap
            0 1 1 1;  # 0 gaps
            1 0 0 1;  # 2 gaps
            0 0 0 0
        ]
        @test FoodWebTools.calculate_gaps(A) == 3

    end

    @testset "intervality" begin

        # Already interval
        A = Bool[
            1 1 0;
            0 1 1;
            0 0 1
        ]

        @test intervality(A) == 0

        # Intervality can never exceed the initial gap count
        A = Bool[
            1 0 1;
            0 1 0;
            1 1 0
        ]

        initial = FoodWebTools.calculate_gaps(A)
        best = intervality(A; iterations=5000)

        @test best <= initial
        @test best >= 0

        # Empty network
        A = falses(5, 5)
        @test intervality(A) == 0

        # Complete network
        A = trues(5, 5)
        @test intervality(A) == 0

    end

end

# downsampling

@testset "Diameter Tests" begin

    @testset "Empty graph" begin
        A = falses(1,1)
        @test diameter(A) == 0
    end

    @testset "Single edge" begin
        # 1 → 2
        A = Bool[
            0 1
            0 0
        ]
        @test diameter(A) == 1
    end

    @testset "Directed chain" begin
        # 1 → 2 → 3 → 4
        A = Bool[
            0 1 0 0
            0 0 1 0
            0 0 0 1
            0 0 0 0
        ]
        @test diameter(A) == 3
    end

    @testset "Directed cycle" begin
        # 1 → 2 → 3 → 1
        A = Bool[
            0 1 0
            0 0 1
            1 0 0
        ]
        @test diameter(A) == 2
    end

    @testset "Disconnected graph" begin
        # 1 → 2     3 → 4
        A = Bool[
            0 1 0 0
            0 0 0 0
            0 0 0 1
            0 0 0 0
        ]
        @test diameter(A) == 1
    end

    @testset "Completely disconnected" begin
        A = falses(5,5)
        @test diameter(A) == 0
    end

    @testset "Complete digraph" begin
        n = 5
        A = trues(n,n)
        A[diagind(A)] .= false

        @test diameter(A) == 1
    end

    @testset "Star graph" begin
        #      1
        #    / | \
        #   2  3  4
        A = Bool[
            0 1 1 1
            0 0 0 0
            0 0 0 0
            0 0 0 0
        ]

        @test diameter(A) == 1
    end

    @testset "Branching graph" begin
        #
        # 1 → 2 → 4
        #  \
        #   → 3 → 5
        #
        A = Bool[
            0 1 1 0 0
            0 0 0 1 0
            0 0 0 0 1
            0 0 0 0 0
            0 0 0 0 0
        ]

        @test diameter(A) == 2
    end

    @testset "Self loops ignored" begin
        A = Bool[
            1 1 0
            0 1 1
            0 0 1
        ]

        @test diameter(A) == 2
    end

end