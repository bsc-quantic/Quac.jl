@testset "Array" begin
    using Quac: I

    @testset "Matrix" begin
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
            Swap,
            Control{I},
            Control{Control{I}},
            Control{Control{Control{I}}},
            Control{Swap},
            Control{Control{Swap}},
            Control{Control{Control{Swap}}},
        ]
            @test Matrix(Gate{Op}) isa Matrix{ComplexF32}
            @test size(Matrix(Gate{Op})) == (2^length(Op), 2^length(Op))
        end

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
            N = length(Op)
            gate = Gate{Op}(1:N...; rand(Op)...)
            @test Matrix(gate) isa Matrix{ComplexF32}
            @test size(Matrix(gate)) == (2^length(Op), 2^length(Op))
        end
    end

    @testset "Diagonal" begin
        using LinearAlgebra: Diagonal, diag

        for Op in [I, Z, S, Sd, T, Td, Control{Z}, Control{Control{Z}}, Control{Control{Control{Z}}}]
            @test Diagonal(Gate{Op}) isa Diagonal{ComplexF32}
        end

        for Op in [
            I,
            Z,
            S,
            Sd,
            T,
            Td,
            Rz,
            Control{Z},
            Control{Control{Z}},
            Control{Control{Control{Z}}},
            Control{Rz},
            Control{Control{Rz}},
            Control{Control{Control{Rz}}},
        ]
            gate = Gate{Op}(1:length(Op)...; rand(Op)...)
            @test Diagonal(gate) isa Diagonal{ComplexF32}
            @test size(Matrix(gate)) == (2^length(Op), 2^length(Op))
        end
    end

    @testset "Array" begin
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
            Swap,
            Control{I},
            Control{Control{I}},
            Control{Control{Control{I}}},
            Control{Swap},
            Control{Control{Swap}},
            Control{Control{Control{Swap}}},
        ]
            N = length(Op) * 2
            @test Array(Gate{Op}) isa Array{ComplexF32,N}
            @test size(Array(Gate{Op})) == tuple(fill(2, N)...)
        end

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
            N = length(Op)
            gate = Gate{Op}(1:N...; rand(Op)...)
            @test Array(gate) isa Array{ComplexF32,2N}
            @test size(Array(gate)) == tuple(fill(2, 2N)...)
        end
    end
end