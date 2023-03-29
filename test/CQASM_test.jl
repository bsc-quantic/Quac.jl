@testset "Serialization" begin
    using Quac.cQASM

    @testset "cQASM_generic" begin

    end

    @testset "cQASM_gates" begin
        @testset "one_qubit" begin
            for entry in [
                "i q[0]" => "i",
                "h q[0]" => "h",
                "x q[0]" => "x",
                "y q[0]" => "y",
                "z q[0]" => "z",
                "x90 q[0]" => "x90",
                "y90 q[0]" => "y90",
                "mx90 q[0]" => "mx90",
                "my90 q[0]" => "my90",
                "s q[0]" => "s",
                "sdag q[0]" => "sdag",
                "t q[0]" => "t",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "i bit1" => [["i", "bit1"]],
                "i quentin2" => [["i", "quentin2"]],
                "i tagname" => [["i", "tagname"]],
                "i _" => [["i", "_"]],
                "i _12" => [["i", "_12"]],
                "i q" => [["i", "q"]],
                "i q4" => [["i", "q4"]],
                "i q[9]" => [["i", "q[9]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "two_qubit" begin
            for entry in [
                "cnot q[0],q[1]" => "cnot",
                "cz q[0],q[1]" => "cz",
                "swap q[0],q[1]" => "swap",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "cnot bit1,bit2" => [["cnot", "bit1", "bit2"]],
                "cnot quentin2,quentin1" => [["cnot", "quentin2", "quentin1"]],
                "cnot tagname,tagname2" => [["cnot", "tagname", "tagname2"]],
                "cnot _,_2" => [["cnot", "_", "_2"]],
                "cnot _12,_122" => [["cnot", "_12", "_122"]],
                "cnot q,q2" => [["cnot", "q", "q2"]],
                "cnot q4,q2" => [["cnot", "q4", "q2"]],
                "cnot q[9],q[2]" => [["cnot", "q[9]", "q[2]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "three_qubit" begin
            for entry in [
                "toffoli q[0],q[1],q[2]" => "toffoli",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "toffoli bit1,bit2,bit3" => [["toffoli", "bit1", "bit2", "bit3"]],
                "toffoli quentin2,quentin1,quentin3" => [["toffoli", "quentin2", "quentin1", "quentin3"]],
                "toffoli tagname,tagname2,tagname3" => [["toffoli", "tagname", "tagname2", "tagname3"]],
                "toffoli _,_2,_3" => [["toffoli", "_", "_2", "_3"]],
                "toffoli _12,_122,_123" => [["toffoli", "_12", "_122", "_123"]],
                "toffoli q,q2,q3" => [["toffoli", "q", "q2", "q3"]],
                "toffoli q4,q2,q3" => [["toffoli", "q4", "q2", "q3"]],
                "toffoli q[9],q[2],q[3]" => [["toffoli", "q[9]", "q[2]", "q[3]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "rotation" begin
            
        end

        @testset "unique" begin
            
        end

        @testset "bit_controlled_one_qubit" begin
            for entry in [
                "c-i b[0],q[1]" => "c-i",
                "c-h b[0],q[1]" => "c-h",
                "c-x b[0],q[1]" => "c-x",
                "c-y b[0],q[1]" => "c-y",
                "c-z b[0],q[1]" => "c-z",
                "c-x90 b[0],q[1]" => "c-x90",
                "c-y90 b[0],q[1]" => "c-y90",
                "c-mx90 b[0],q[1]" => "c-mx90",
                "c-my90 b[0],q[1]" => "c-my90",
                "c-s b[0],q[1]" => "c-s",
                "c-sdag b[0],q[1]" => "c-sdag",
                "c-t b[0],q[1]" => "c-t",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "c-i bit1,q[1]" => [["c-i", "bit1", "q[1]"]],
                "c-i quentin2,q[1]" => [["c-i", "quentin2", "q[1]"]],
                "c-i tagname,q[1]" => [["c-i", "tagname", "q[1]"]],
                "c-i _,q[1]" => [["c-i", "_", "q[1]"]],
                "c-i _12,q[1]" => [["c-i", "_12", "q[1]"]],
                "c-i q,q[1]" => [["c-i", "q", "q[1]"]],
                "c-i q4,q[1]" => [["c-i", "q4", "q[1]"]],
                "c-i b[0],q[1]" => [["c-i", "b[0]", "q[1]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "bit_controlled_two_qubit" begin
            
        end

        @testset "bit_controlled_three_qubit" begin
            
        end

        @testset "bit_controlled_rotation" begin
            
        end

        @testset "bit_controlled_unique" begin
            
        end

    end

    @testset "cQASM_misc" begin
        
    end
end