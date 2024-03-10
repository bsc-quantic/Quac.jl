@testset "Operator" begin
    using Quac: I

    @testset "isparametric" begin
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Swap, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
            @test isparametric(Op) == false
        end

        for Op in [Rx, Ry, Rz, U2, U3, Control{Rx}, Control{Control{Rx}}, Control{Control{Control{Rx}}}]
            @test isparametric(Op) == true
        end
    end

    @testset "parameters" begin
        @test parameters(Control{U2}) === parameters(U2)
        @test parameters(Control{Control{U2}}) === parameters(U2)
        @test parameters(Control{Control{Control{U2}}}) === parameters(U2)
    end

    @testset "Base.length" begin
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz]
            @test length(Op) == 1
        end

        @test length(Swap) == 2

        @test length(Control{I}) == 2
        @test length(Control{Control{I}}) == 3
        @test length(Control{Control{Control{I}}}) == 4
    end

    @testset "Base.adjoint" begin
        for Op in [I, X, Y, Z, H, Rx, Ry, Rz, Swap, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
            @test adjoint(Op()) isa Op
        end

        @test adjoint(S()) == Sd()
        @test adjoint(Sd()) == S()
        @test adjoint(T()) == Td()
        @test adjoint(Td()) == T()
        @test adjoint(Control{S}()) === Control{Sd}()
        @test adjoint(Control{Control{S}}()) === Control{Control{Sd}}()
        @test adjoint(Control{Control{Control{S}}}()) === Control{Control{Control{Sd}}}()
    end

    @testset "Base.rand" begin
        for Op in [
            I,
            X,
            Y,
            Z,
            H,
            S,
            Sd,
            T,
            Td,
            Rx,
            Ry,
            Rz,
            Swap,
            Control{I},
            Control{Control{I}},
            Control{Control{Control{I}}},
            Control{Rx},
            Control{Control{Rx}},
            Control{Control{Control{Rx}}},
            Control{Swap},
            Control{Control{Swap}},
            Control{Control{Control{Swap}}},
        ]
            @test rand(Op) isa Op
        end
    end

    @testset "Operator sets" begin
        @testset "Pauli" begin
            for Op in [I, X, Y, Z]
                @test Op <: Pauli
            end

            for Op in [Rx, Ry, Rz, S, Sd, T, Td, H, Swap, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
                @test !(Op <: Pauli)
            end
        end

        @testset "Phase" begin
            for Op in [I, Z, S, Sd, T, Td, Rz]
                @test Op <: Phase
            end

            for Op in [X, Y, H, Rx, Ry, Swap, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
                @test !(Op <: Phase)
            end
        end
    end
end
