@testset "Array" begin
    using Quac: I
    import LinearAlgebra

    @testset "Matrix" begin
        @testset "$Op()" for Op in [
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
            N = length(Op)
            @test Matrix(Op()) isa Matrix{ComplexF32}
            @test size(Matrix(Op())) == (2^N, 2^N)
        end

        @testset "rand($Op)" for Op in [
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
            Rxx,
            Ryy,
            Rzz,
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
            op = rand(Op)
            @test Matrix(op) isa Matrix{ComplexF32}
            @test size(Matrix(op)) == (2^length(Op), 2^length(Op))
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            mat = Matrix(rand(SU{N}))
            @test mat isa Matrix{ComplexF64}
            @test size(mat) == (2^N, 2^N)
        end
    end

    @testset "Diagonal" begin
        using LinearAlgebra: Diagonal

        @testset "$Op" for Op in [
            I,
            Z,
            S,
            Sd,
            T,
            Td,
            Rz,
            Rzz,
            Control{Z},
            Control{Control{Z}},
            Control{Control{Control{Z}}},
            Control{Rz},
            Control{Control{Rz}},
            Control{Control{Control{Rz}}},
        ]
            N = length(Op)
            op = rand(Op)
            @test Diagonal(op) isa Diagonal{ComplexF32}
            @test size(Diagonal(op)) == (2^N, 2^N)
        end
    end

    @testset "Array" begin
        @testset "$Op()" for Op in [
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
            N = length(Op)
            op = rand(Op)
            @test Array(op) isa Array{ComplexF32,2N}
            @test size(Array(op)) == tuple(fill(2, 2N)...)
        end

        @testset "rand($Op)" for Op in [
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
            Rxx,
            Ryy,
            Rzz,
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
            op = rand(Op)
            @test Array(gate) isa Array{ComplexF32,2N}
            @test size(Array(gate)) == tuple(fill(2, 2N)...)
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            op = rand(SU{N})
            @test Array(op) isa Array{ComplexF32,2N}
            @test size(Array(op)) == tuple(fill(2, 2N)...)
        end
    end

    @testset "Parametric operators" begin
        @testset "Rxx" begin
            @test begin
                op = Rxx(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Rxx(; θ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 0 0 -1im; 0 1 -1im 0; 0 -1im 1 0; -1im 0 0 1]
            end

            @test begin
                op = Rxx(; θ = π)
                x = Matrix(X())
                Matrix(op) ≈ -1im * reshape(permutedims(reshape(kron(vec(x), vec(x)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end

        @testset "Ryy" begin
            @test begin
                op = Ryy(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Ryy(; θ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 0 0 1im; 0 1 -1im 0; 0 -1im 1 0; 1im 0 0 1]
            end

            @test begin
                op = Ryy(; θ = π)
                y = Matrix(Y())
                Matrix(op) ≈ -1im * reshape(permutedims(reshape(kron(vec(y), vec(y)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end

        @testset "Rzz" begin
            @test begin
                op = Rzz(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Rzz(; θ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * Diagonal([1 - 1im, 1 + 1im, 1 + 1im, 1 - 1im])
            end

            @test begin
                op = Rzz(; θ = π)
                z = Matrix(Z())
                Matrix(op) ≈ -1im * reshape(permutedims(reshape(kron(vec(z), vec(z)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end
    end
end
