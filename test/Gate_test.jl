@testset "Gate" begin
    using Quac: I
    using LinearAlgebra: qr

    @testset "Base.length" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
            Control{I},
            Control{Control{I}},
            Control{Control{Control{I}}},
            Control{Swap},
            Control{Control{Swap}},
            Control{Control{Control{Swap}}},
            SU{2},
            SU{4},
        ]
            @test length(Gate{Op}) == length(Op)
        end

        @testset "Product operator (⊗)" begin
            @test length(X() ⊗ Y()) == 2
            @test length(X() ⊗ Y() ⊗ Z()) == 3
            @test length(X() ⊗ Y() ⊗ Z() ⊗ H()) == 4
            @test length(X() ⊗ Swap()) == 3
            @test length(Swap() ⊗ CRx()) == 4
        end
    end

    @testset "Constructor" begin
        @testset "$Op" for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3]
            N = length(Op)
            @test Gate{Op}(1:N...; Quac.randtuple(parameters(Op))...) |> !isnothing
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            q, _ = qr(rand(ComplexF32, 2^N, 2^N))
            matrix = Matrix{ComplexF32}(q)
            @test Gate{SU{N}}(1:N...; matrix) |> !isnothing
        end

        @testset "SU - unitary condition" begin
            # test_throws on a non-unitary matrix
            @test_throws ArgumentError SU{2}(1, 2; matrix = rand(ComplexF32, 4, 4))

            # test_throws on a non-square matrix
            @test_throws ArgumentError SU{2}(1, 2; matrix = rand(ComplexF32, 4, 2))

            # test_throws on a matrix without size (N, N)
            @test_throws ArgumentError SU{2}(1, 2; matrix = rand(ComplexF32, 2, 2))

            # test_throws SU{N} with N not a power of 2
            @test_throws ArgumentError SU{2}(1, 2; matrix = rand(ComplexF32, 3, 3))

            # test_throws when there are not N lanes
            q, _ = qr(rand(ComplexF32, 4, 4))
            @test_throws ArgumentError SU{3}(1, 2, 3; matrix = Matrix{ComplexF32}(q))
        end

        @testset "Product operator (⊗)" begin
            op = X() ⊗ Y()
            N = length(op)
            @test op(1:N...) isa Gate{Quac.ProductOperator,N}

            op = X() ⊗ Y() ⊗ Z()
            N = length(op)
            @test op(1:N...) isa Gate{Quac.ProductOperator,N}

            op = X() ⊗ Y() ⊗ Z() ⊗ H()
            N = length(op)
            @test op(1:N...) isa Gate{Quac.ProductOperator,N}

            op = X() ⊗ Swap()
            N = length(op)
            @test op(1:N...) isa Gate{Quac.ProductOperator,N}

            op = Swap() ⊗ CRx()
            N = length(op)
            @test op(1:N...) isa Gate{Quac.ProductOperator,N}
        end
    end

    @testset "Constructor aliases" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            @test Op(1:N...; Quac.randtuple(parameters(Op))...) isa Gate{Op,N}
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            q, _ = qr(rand(ComplexF64, 2^N, 2^N))
            matrix = Matrix{ComplexF64}(q)
            @test SU{N}(1:N...; matrix) isa Gate{SU{N},N}
        end
    end

    @testset "lanes" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            @test lanes(Gate{Op}(1:N...; Quac.randtuple(parameters(Op))...)) === tuple(1:N...)
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            q, _ = qr(rand(ComplexF32, 2^N, 2^N))
            matrix = Matrix{ComplexF32}(q)
            @test lanes(Gate{SU{N}}(1:N...; matrix)) === tuple(1:length(SU{N})...)
        end
    end

    @testset "operator" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            # NOTE `operator(Gate{Op})` is implied
            @test operator(Gate{Op}(1:length(Op)...; Quac.randtuple(parameters(Op))...)) isa Op
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            q, _ = qr(rand(ComplexF32, 2^N, 2^N))
            matrix = Matrix{ComplexF32}(q)
            @test operator(Gate{SU{N}}(1:N...; matrix)) isa SU{N}
        end
    end

    @testset "parameters" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            @test parameters(Gate{Op}) === parameters(Op)

            params = Quac.randtuple(parameters(Op))
            @test parameters(Gate{Op}(1:length(Op)...; params...)) == params
        end

        # Special case for SU{N}
        @testset "SU{$N}" for N in [1, 2, 3]
            @test parameters(Gate{SU{N}}) === parameters(SU{N})

            q, _ = qr(rand(ComplexF32, 2^N, 2^N))
            matrix = Matrix{ComplexF32}(q)
            @test parameters(Gate{SU{N}}(1:N...; matrix)).matrix == matrix
        end
    end

    @testset "adjoint" begin
        @testset "Gate type - $Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            @test adjoint(Gate{Op}) === Gate{adjoint(Op),length(Op)}
        end

        # `adjoint(::Gate)` with no parameters
        @testset "Gate instance - $Op" for Op in [
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
        ]
            N = length(Op)
            @test adjoint(Op(1:N...)) === Op'(1:N...)
            @test adjoint(Gate{Op}(1:N...)) === Gate{Op',N}(1:N...)
        end

        # `adjoint(::Gate)` with parameters
        @testset "Gate instance - $Op" for Op in [
            Rx,
            Ry,
            Rz,
            Rxx,
            Ryy,
            Rzz,
            U2,
            U3,
            Control{Rx},
            Control{Control{Rx}},
            Control{Control{Control{Rx}}},
        ]
            params = Quac.randtuple(parameters(Op))
            N = length(Op)
            @test Op(1:N...; params...)' === Op'(1:N...; [key => -val for (key, val) in pairs(params)]...)
            @test adjoint(Gate{Op}(1:N...; params...)) ===
                  Gate{adjoint(Op),N}(1:N...; [key => -val for (key, val) in pairs(params)]...)
        end

        # Special case for SU{N}
        @testset "Gate instance - SU{$N}" for N in [1, 2, 3]
            q, _ = qr(rand(ComplexF64, 2^N, 2^N))
            q = Matrix(q)

            @test adjoint(SU{N}(1:N...; matrix = q)) == SU{N}'(1:N...; matrix = q')
            @test adjoint(Gate{SU{N}}(1:N...; matrix = q)) == Gate{SU{N}'}(1:N...; matrix = q')
        end
    end

    @testset "Base.rand" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
            Control{I},
            Control{Control{I}},
            Control{Control{Control{I}}},
            Control{Rx},
            Control{Control{Rx}},
            Control{Control{Control{Rx}}},
            Control{Swap},
            Control{Control{Swap}},
            Control{Control{Control{Swap}}},
            SU{2},
            SU{4},
            SU{8},
        ]
            N = length(Op)
            @test rand(Gate{Op}, 1:N...) isa Gate{Op,N}
        end
    end

    @testset "Base.propertynames" begin
        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
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
            @test Quac.ntnames(parameters(Op)) ⊆ propertynames(Gate{Op})

            N = length(Op)
            @test Quac.ntnames(parameters(Op)) ⊆ propertynames(rand(Gate{Op}, 1:N...))
        end
    end

    @testset "Base.getproperty" begin
        # TODO
    end

    @testset "targettype" begin
        using Quac: targettype

        @testset "$Op" for Op in [
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
            U2,
            U3,
            Swap,
            FSim,
            SU{2},
            SU{4},
            SU{8},
        ]
            @test targettype(Gate{Op}) === Op
        end

        @testset "$Op" for Op in [I, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
            @test targettype(Gate{Op}) === I
        end

        @testset "$(operator(gate))" for (gate, correct) in [
            (Z(1), Z),
            (CZ(1, 2), Z),
            (Control(1, Control(2, Z(3))), Z),
        ]
            @test targettype(gate) === correct
            @test targettype(operator(gate)) === correct
        end
    end

    @testset "control" begin
        @testset "throw error - $Op" for Op in
                                         [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Swap, FSim]
            @test_throws MethodError control(rand(Gate{Op}, 1:length(Op)...))
        end

        @testset "$Op" for Op in [I, Rx, Swap, Rxx, U2, U3, FSim]
            N = length(Op)

            M = length(Control{Op})
            @test control(rand(Gate{Control{Op}}, 1:M...)) === (1:(M-N)...,)

            M = length(Control{Control{Op}})
            @test control(rand(Gate{Control{Control{Op}}}, 1:M...)) === (1:(M-N)...,)

            M = length(Control{Control{Control{Op}}})
            @test control(rand(Gate{Control{Control{Control{Op}}}}, 1:M...)) === (1:(M-N)...,)
        end
    end

    @testset "target" begin
        @testset "throw error - $Op" for Op in
                                         [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Swap, FSim]
            @test_throws MethodError target(rand(Gate{Op}, 1:length(Op)...))
        end

        @testset "$Op" for Op in [I, Rx, Swap, Rxx, U2, U3, FSim]
            N = length(Op)

            M = length(Control{Op})
            @test target(rand(Gate{Control{Op}}, 1:M...)) === ((M-N)+1:M...,)

            M = length(Control{Control{Op}})
            @test target(rand(Gate{Control{Control{Op}}}, 1:M...)) === ((M-N)+1:M...,)

            M = length(Control{Control{Control{Op}}})
            @test target(rand(Gate{Control{Control{Control{Op}}}}, 1:M...)) === ((M-N)+1:M...,)
        end
    end
end
