@testset "Gate" begin
    using Quac: I

    @testset "Base.length" begin
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
        ]
            @test length(Gate{Op}) === length(Op)
        end
    end

    @testset "Constructor" begin
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3]
            @test Gate{Op}(range(1, length = length(Op))...; rand(parameters(Op))...) !== nothing
        end
    end

    @testset "Constructor aliases" begin
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
            P = parameters(Op)
            @test Op(1:N...; rand(P)...) isa Gate{Op,N,P}
        end
    end

    @testset "lanes" begin
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
            @test lanes(Gate{Op}(1:length(Op)...; rand(parameters(Op))...)) === tuple(1:length(Op)...)
        end
    end

    @testset "operator" begin
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
            @test operator(Gate{Op}(1:length(Op)...; rand(parameters(Op))...)) === Op
        end
    end

    @testset "parameters" begin
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

            params = rand(parameters(Op))
            @test parameters(Gate{Op}(1:length(Op)...; params...)) === params
        end
    end

    @testset "adjoint" begin
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
            @test adjoint(Gate{Op}) === Gate{adjoint(Op)}
        end

        # `adjoint(::Gate)` with no parameters
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Swap, Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
            @test adjoint(Gate{Op}(1:length(Op)...)) === Gate{adjoint(Op)}(1:length(Op)...)
        end

        # `adjoint(::Gate)` with parameters
        for Op in [Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Control{Rx}, Control{Control{Rx}}, Control{Control{Control{Rx}}}]
            params = rand(parameters(Op))
            @test adjoint(Gate{Op}(1:length(Op)...; params...)) ===
                  Gate{adjoint(Op)}(1:length(Op)...; [key => -val for (key, val) in pairs(params)]...)
        end
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
            P = parameters(Op)
            @test rand(Gate{Op}, 1:N...) isa Gate{Op,N,P}
        end
    end

    @testset "Base.propertynames" begin
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
            @test keys(parameters(Op)) ⊆ propertynames(Gate{Op})

            N = length(Op)
            @test keys(parameters(Op)) ⊆ propertynames(rand(Gate{Op}, 1:N...))
        end
    end

    @testset "Base.getproperty" begin
        # TODO
    end

    @testset "targettype" begin
        using Quac: targettype

        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Swap, FSim]
            @test targettype(Gate{Op}) === Op
        end

        for Op in [Control{I}, Control{Control{I}}, Control{Control{Control{I}}}]
            @test targettype(Gate{Op}) === I
        end
    end

    @testset "control" begin
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Swap, FSim]
            @test_throws MethodError control(rand(Gate{Op}, 1:length(Op)...))
        end

        for Op in [I, Rx, Swap, Rxx, U2, U3, FSim]
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
        for Op in [I, X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Rxx, Ryy, Rzz, U2, U3, Swap, FSim]
            @test_throws MethodError target(rand(Gate{Op}, 1:length(Op)...))
        end

        for Op in [I, Rx, Swap, Rxx, U2, U3, FSim]
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
