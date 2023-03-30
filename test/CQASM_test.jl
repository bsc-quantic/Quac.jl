@testset "Serialization" begin
    using Quac.cQASM

    @testset "cQASM_generic" begin
        # for entry in [
        #     "version 1" => [["version", "1"]],
        #     "version 2.3" => [["version", "2.3"]],
        #     "qubits 4" => [["qubits", "4"]],
        #     "qubits 52" => [["qubits", "52"]],
        #     "i tagname" => [["i", "tagname"]],
        #     "i _" => [["i", "_"]],
        #     "i _12" => [["i", "_12"]],
        # ]
        #     @test parseCQASMCode(entry.first) == entry.second
        # end
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
            for entry in [
                "rx q[0],3.14" => "rx",
                "ry q[0],3.14" => "ry",
                "rz q[0],3.14" => "rz",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "rx bit1,0.71" => [["rx", "bit1", "0.71"]],
                "rx quentin2,0.71" => [["rx", "quentin2", "0.71"]],
                "rx tagname,0.71" => [["rx", "tagname", "0.71"]],
                "rx _,0.71" => [["rx", "_", "0.71"]],
                "rx _12,0.71" => [["rx", "_12", "0.71"]],
                "rx q,0.71" => [["rx", "q", "0.71"]],
                "rx q4,0.71" => [["rx", "q4", "0.71"]],
                "rx q[9],0.71" => [["rx", "q[9]", "0.71"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "unique" begin
            for entry in [
                "crk q[0],q[1],5" => "crk",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "crk bit1,bit2,3" => [["crk", "bit1", "bit2", "3"]],
                "crk quentin2,quentin1,3" => [["crk", "quentin2", "quentin1", "3"]],
                "crk tagname,tagname2,3" => [["crk", "tagname", "tagname2", "3"]],
                "crk _,_2,3" => [["crk", "_", "_2", "3"]],
                "crk _12,_122,3" => [["crk", "_12", "_122", "3"]],
                "crk q,q2,3" => [["crk", "q", "q2", "3"]],
                "crk q4,q2,3" => [["crk", "q4", "q2", "3"]],
                "crk q[9],q[2],3" => [["crk", "q[9]", "q[2]", "3"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
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
            for entry in [
                "c-cnot b[3],q[0],q[1]" => "c-cnot",
                "c-cz b[3],q[0],q[1]" => "c-cz",
                "c-swap b[3],q[0],q[1]" => "c-swap",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "c-cnot bit0,bit1,bit2" => [["c-cnot", "bit0", "bit1", "bit2"]],
                "c-cnot quentin0,quentin2,quentin1" => [["c-cnot", "quentin0", "quentin2", "quentin1"]],
                "c-cnot tagname0,tagname,tagname2" => [["c-cnot", "tagname0", "tagname", "tagname2"]],
                "c-cnot _,_0,_2" => [["c-cnot", "_", "_0", "_2"]],
                "c-cnot _120,_12,_122" => [["c-cnot", "_120", "_12", "_122"]],
                "c-cnot q,q0,q2" => [["c-cnot", "q", "q0", "q2"]],
                "c-cnot q40,q4,q2" => [["c-cnot", "q40", "q4", "q2"]],
                "c-cnot b[0],q[9],q[2]" => [["c-cnot", "b[0]", "q[9]", "q[2]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "bit_controlled_three_qubit" begin
            for entry in [
                "c-toffoli b[4],q[0],q[1],q[2]" => "c-toffoli",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "c-toffoli bit0,bit1,bit2,bit3" => [["c-toffoli", "bit0", "bit1", "bit2", "bit3"]],
                "c-toffoli quentin0,quentin2,quentin1,quentin3" => [["c-toffoli", "quentin0", "quentin2", "quentin1", "quentin3"]],
                "c-toffoli tagname0,tagname,tagname2,tagname3" => [["c-toffoli", "tagname0", "tagname", "tagname2", "tagname3"]],
                "c-toffoli _,_0,_2,_3" => [["c-toffoli", "_", "_0", "_2", "_3"]],
                "c-toffoli _120,_12,_122,_123" => [["c-toffoli", "_120", "_12", "_122", "_123"]],
                "c-toffoli q,q0,q2,q3" => [["c-toffoli", "q", "q0", "q2", "q3"]],
                "c-toffoli q40,q4,q2,q3" => [["c-toffoli", "q40", "q4", "q2", "q3"]],
                "c-toffoli b[0],q[9],q[2],q[3]" => [["c-toffoli", "b[0]", "q[9]", "q[2]", "q[3]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "bit_controlled_rotation" begin
            for entry in [
                "c-rx b[1],q[0],3.14" => "c-rx",
                "c-ry b[1],q[0],3.14" => "c-ry",
                "c-rz b[1],q[0],3.14" => "c-rz",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "c-rx bit0,bit1,0.71" => [["c-rx", "bit0", "bit1", "0.71"]],
                "c-rx quentin0,quentin2,0.71" => [["c-rx", "quentin0", "quentin2", "0.71"]],
                "c-rx tagname0,tagname,0.71" => [["c-rx", "tagname0", "tagname", "0.71"]],
                "c-rx _,_0,0.71" => [["c-rx", "_", "_0", "0.71"]],
                "c-rx _12,_120,0.71" => [["c-rx", "_12", "_120", "0.71"]],
                "c-rx q,q0,0.71" => [["c-rx", "q", "q0", "0.71"]],
                "c-rx q4,q40,0.71" => [["c-rx", "q4", "q40", "0.71"]],
                "c-rx b[0],q[9],0.71" => [["c-rx", "b[0]", "q[9]", "0.71"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "bit_controlled_unique" begin
            for entry in [
                "c-crk b[0],q[0],q[1],5" => "c-crk",
            ]
                @test parseCQASMCode(entry.first)[1][1] == entry.second
            end

            for entry in [
                "c-crk bit0,bit1,bit2,3" => [["c-crk", "bit0", "bit1", "bit2", "3"]],
                "c-crk quentin0,quentin2,quentin1,3" => [["c-crk", "quentin0", "quentin2", "quentin1", "3"]],
                "c-crk tagname0,tagname,tagname2,3" => [["c-crk", "tagname0", "tagname", "tagname2", "3"]],
                "c-crk _,_0,_2,3" => [["c-crk", "_", "_0", "_2", "3"]],
                "c-crk _12,_120,_122,3" => [["c-crk", "_12", "_120", "_122", "3"]],
                "c-crk q,q0,q2,3" => [["c-crk", "q", "q0", "q2", "3"]],
                "c-crk q4,q40,q2,3" => [["c-crk", "q4", "q40", "q2", "3"]],
                "c-crk b[0],q[9],q[2],3" => [["c-crk", "b[0]", "q[9]", "q[2]", "3"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

    end

    @testset "cQASM_misc" begin
        
    end
end