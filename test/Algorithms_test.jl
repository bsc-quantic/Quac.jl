@testset "Algorithms" begin

    @testset "Quac.Algorithms.QFT" begin
        circuit = Quac.Algorithms.QFT(3)

        @test length(circuit.lanes) == 3
    end

    @testset "Quac.Algorithms.QuantumVolume" begin
        n_qubits = 4
        depth = 2

        circuit = Quac.Algorithms.QuantumVolume(n_qubits, depth)

        @test length(circuit.lanes) == n_qubits
        @test length(circuit) == (n_qubits - n_qubits % 2) * depth
    end
end