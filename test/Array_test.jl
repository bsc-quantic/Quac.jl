@testset "Array" begin
    using Quac: I
    import LinearAlgebra

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

        # Special case for SU{N}
        for N in [2, 4, 8]
            @test Matrix(rand(Gate{SU{N}}, 1:Int(log2(N))...)) isa Matrix{ComplexF32}
            @test size(Matrix(rand(Gate{SU{N}}, 1:Int(log2(N))...))) == (N, N)
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
            Rzz,
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

        # Special case for SU{N}
        for N in [2, 4, 8]
            n_qubits = Int(log2(N))
            @test Array(rand(Gate{SU{N}}, 1:Int(log2(N))...)) isa Array{ComplexF32, 2*n_qubits}
            @test size(Array(rand(Gate{SU{N}}, 1:Int(log2(N))...))) == tuple(fill(2, 2*n_qubits)...)
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
            gate = Gate{Op}(1:N...; rand(Op)...)
            @test Array(gate) isa Array{ComplexF32,2N}
            @test size(Array(gate)) == tuple(fill(2, 2N)...)
        end
    end

    @testset "Correctness" begin
        @testset "Rxx" begin
            @test begin
                g = Rxx(1, 2, θ = 0)
                Matrix(g) ≈ LinearAlgebra.I(2^length(g))
            end

            @test begin
                g = Rxx(1, 2, θ = π / 2)
                Matrix(g) ≈ 1 / sqrt(2) * [1 0 0 -1im; 0 1 -1im 0; 0 -1im 1 0; -1im 0 0 1]
            end

            @test begin
                g = Rxx(1, 2, θ = π)
                x = Matrix(Gate{X})
                Matrix(g) ≈ -1im * reshape(permutedims(reshape(kron(vec(x), vec(x)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end

        @testset "Ryy" begin
            @test begin
                g = Ryy(1, 2, θ = 0)
                Matrix(g) ≈ LinearAlgebra.I(2^length(g))
            end

            @test begin
                g = Ryy(1, 2, θ = π / 2)
                Matrix(g) ≈ 1 / sqrt(2) * [1 0 0 1im; 0 1 -1im 0; 0 -1im 1 0; 1im 0 0 1]
            end

            @test begin
                g = Ryy(1, 2, θ = π)
                x = Matrix(Gate{Y})
                Matrix(g) ≈ -1im * reshape(permutedims(reshape(kron(vec(x), vec(x)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end

        @testset "Rzz" begin
            @test begin
                g = Rzz(1, 2, θ = 0)
                Matrix(g) ≈ LinearAlgebra.I(2^length(g))
            end

            @test begin
                g = Rzz(1, 2, θ = π / 2)
                Matrix(g) ≈ 1 / sqrt(2) * Diagonal([1 - 1im, 1 + 1im, 1 + 1im, 1 - 1im])
            end

            @test begin
                g = Rzz(1, 2, θ = π)
                x = Matrix(Gate{Z})
                Matrix(g) ≈ -1im * reshape(permutedims(reshape(kron(vec(x), vec(x)), 2, 2, 2, 2), (1, 3, 2, 4)), 4, 4)
            end
        end
    end
end
