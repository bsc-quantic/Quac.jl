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
            @test Matrix(Op()) isa Matrix{ComplexF64}
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
            @test Matrix(op) isa Matrix{ComplexF64}
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
            @test Diagonal(op) isa Diagonal{ComplexF64}
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
            @test Array(op) isa Array{ComplexF64,2N}
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
            SU{1},
            SU{2},
            SU{3},
        ]
            N = length(Op)
            op = rand(Op)
            @test Array(op) isa Array{ComplexF64,2N}
            @test size(Array(op)) == tuple(fill(2, 2N)...)
        end
    end

    @testset "Parametric operators" begin
        @testset "Rx" begin
            @test begin
                op = Rx(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Rx(; θ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1im; -1im 1]
            end

            @test begin
                op = Rx(; θ = π)
                Matrix(op) ≈ [0 -1im; -1im 0]
            end
        end

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

        @testset "Ry" begin
            @test begin
                op = Ry(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Ry(; θ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1; 1 1]
            end

            @test begin
                op = Ry(; θ = π)
                Matrix(op) ≈ [0 -1; 1 0]
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
                Matrix(op) ≈ -1im * kron(y, y)
            end
        end

        @testset "Rz" begin
            @test begin
                op = Rz(; θ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Rz(; θ = π / 2)
                Matrix(op) ≈ [1 0; 0 1im]
            end

            @test begin
                op = Rz(; θ = π)
                Matrix(op) ≈ [1 0; 0 -1]
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

        @testset "U2" begin
            @test begin
                op = U2(; ϕ = 0, λ = 0)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1; 1 1]
            end

            @test begin
                op = U2(; ϕ = 0, λ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1im; 1 1im]
            end

            @test begin
                op = U2(; ϕ = 0, λ = π)
                Matrix(op) ≈ 1 / sqrt(2) * [1 1; 1 -1]
            end

            @test begin
                op = U2(; ϕ = π / 2, λ = 0)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1; 1im 1im]
            end

            @test begin
                op = U2(; ϕ = π / 2, λ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1im; 1im -1]
            end

            @test begin
                op = U2(; ϕ = π / 2, λ = π)
                Matrix(op) ≈ 1 / sqrt(2) * [1 1; 1im -1im]
            end

            @test begin
                op = U2(; ϕ = π, λ = 0)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1; -1 -1]
            end

            @test begin
                op = U2(; ϕ = π, λ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1im; -1 -1im]
            end

            @test begin
                op = U2(; ϕ = π, λ = π)
                Matrix(op) ≈ 1 / sqrt(2) * [1 1; -1 1]
            end
        end

        @testset "U3" begin
            @test begin
                op = U3(; θ = 0, ϕ = 0, λ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = U3(; θ = π / 2, ϕ = π / 2, λ = π / 2)
                Matrix(op) ≈ 1 / sqrt(2) * [1 -1im; 1im -1]
            end

            @test begin
                op = U3(; θ = π, ϕ = π, λ = π)
                Matrix(op) ≈ [0 1; -1 0]
            end
        end

        @testset "Hz" begin
            @test begin
                op = Hz(; θ = 0, ϕ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = Hz(; θ = π / 2, ϕ = π / 2)
                Matrix(op) ≈ 1 / 2 * [1+1im 1-1im; -1+1im 1+1im]
            end

            @test begin
                op = Hz(; θ = π, ϕ = π)
                Matrix(op) ≈ [0 1; 1 0]
            end
        end

        @testset "FSim" begin
            @test begin
                op = FSim(; θ = 0, ϕ = 0)
                N = length(op)
                Matrix(op) ≈ LinearAlgebra.I(2^N)
            end

            @test begin
                op = FSim(; θ = 0, ϕ = π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1im]
            end

            @test begin
                op = FSim(; θ = 0, ϕ = π)
                Matrix(op) ≈ [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]
            end

            @test begin
                op = FSim(; θ = 0, ϕ = 3 * π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1im]
            end

            @test begin
                op = FSim(; θ = π / 2, ϕ = 0)
                Matrix(op) ≈ [1 0 0 0; 0 0 -1im 0; 0 -1im 0 0; 0 0 0 1]
            end

            @test begin
                op = FSim(; θ = π / 2, ϕ = π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 0 -1im 0; 0 -1im 0 0; 0 0 0 -1im]
            end

            @test begin
                op = FSim(; θ = π / 2, ϕ = π)
                Matrix(op) ≈ [1 0 0 0; 0 0 -1im 0; 0 -1im 0 0; 0 0 0 -1]
            end

            @test begin
                op = FSim(; θ = π / 2, ϕ = 3 * π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 0 -1im 0; 0 -1im 0 0; 0 0 0 1im]
            end

            @test begin
                op = FSim(; θ = π, ϕ = 0)
                Matrix(op) ≈ [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
            end

            @test begin
                op = FSim(; θ = π, ϕ = π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1im]
            end

            @test begin
                op = FSim(; θ = π, ϕ = π)
                Matrix(op) ≈ [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]
            end

            @test begin
                op = FSim(; θ = π, ϕ = 3 * π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1im]
            end

            @test begin
                op = FSim(; θ = 3 * π / 2, ϕ = 0)
                Matrix(op) ≈ [1 0 0 0; 0 0 1im 0; 0 1im 0 0; 0 0 0 1]
            end

            @test begin
                op = FSim(; θ = 3 * π / 2, ϕ = π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 0 1im 0; 0 1im 0 0; 0 0 0 -1im]
            end

            @test begin
                op = FSim(; θ = 3 * π / 2, ϕ = π)
                Matrix(op) ≈ [1 0 0 0; 0 0 1im 0; 0 1im 0 0; 0 0 0 -1]
            end

            @test begin
                op = FSim(; θ = 3 * π / 2, ϕ = 3 * π / 2)
                Matrix(op) ≈ [1 0 0 0; 0 0 1im 0; 0 1im 0 0; 0 0 0 1im]
            end
        end
    end
end
