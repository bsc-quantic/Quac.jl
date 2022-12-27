using Quac: Circuit, lanes, X

@testset "Circuit" verbose = true begin
    @testset "lanes" begin
        for n in [0, 1, 2, 4]
            @test lanes(Circuit(n)) == n
        end
    end

    @testset "Base.length" begin
        for n in [0, 1, 2]
            @test length(Circuit(n)) == 0
        end
    end

    @testset "Base.isempty" begin
        @test isempty(Circuit(0))

        c = Circuit(1)
        @test isempty(c)

        push!(c, X(1))
        @test !isempty(c)
    end

    @testset "Base.adjoint" begin end

    @testset "Base.vcat" begin end

    @testset "Base.hcat" begin end
end