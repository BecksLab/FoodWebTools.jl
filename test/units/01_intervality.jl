using FoodWebTools
using Test

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
        A = falses(4,4)
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
        A = falses(5,5)
        @test intervality(A) == 0

        # Complete network
        A = trues(5,5)
        @test intervality(A) == 0

    end

end