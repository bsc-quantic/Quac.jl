@testset "Serialization" begin
    using Quac.cQASM

    @testset "cQASM_generic" begin
        for entry in [
            "version 1" => [["version", "1"]],
            "version 2.3" => [["version", "2.3"]],
            "qubits 4" => [["qubits", "4"]],
            "qubits 52" => [["qubits", "52"]],
            "map q[0],data" => [["map", "q[0]", "data"]],
            "map b[1],bit1" => [["map", "b[1]", "bit1"]],
            "prep_z q[3]" => [["prep_z", "q[3]"]],
            "prep_x q[12]" => [["prep_x", "q[12]"]],
            "measure_x q[4]" => [["measure_x", "q[4]"]],
            "measure_z q[32]" => [["measure_z", "q[32]"]],
            "measure q[5]" => [["measure", "q[5]"]],
            "measure_all" => [["measure_all"]],
            "measure_parity q[6],y,q[7],z" => [["measure_parity", "q[6]", "y", "q[7]", "z"]],
            "display" => [["display"]],
            "display b[0]" => [["display", "b[0]"]],
            ".init" => [[".init"]],
            ".grover(3)" => [[".grover", "3"]]
        ]
            @test parseCQASMCode(entry.first) == entry.second
        end

        @test parseCQASMCode("# this is a comment") == Vector{String}[]
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
                "cnot qubit9,q[2]" => [["cnot", "qubit9", "q[2]"]],
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
                "c-i b,q[1]" => [["c-i", "b", "q[1]"]],
                "c-i b1,q[1]" => [["c-i", "b1", "q[1]"]],
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
                "c-cnot b,b,bit2" => [["c-cnot", "b", "b", "bit2"]],
                "c-cnot b,b1,bit2" => [["c-cnot", "b", "b1", "bit2"]],
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
                "c-toffoli b,b,b2,bit3" => [["c-toffoli", "b", "b", "b2", "bit3"]],
                "c-toffoli bit0,bit1,bit2,bit3" => [["c-toffoli", "bit0", "bit1", "bit2", "bit3"]],
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
                "c-rx b,b1,0.71" => [["c-rx", "b", "b1", "0.71"]],
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
                "c-crk b,b1,b,3" => [["c-crk", "b", "b1", "b", "3"]],
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

    @testset "cQASM_extreme_cases" begin
        @testset "multi_binary_controlled" begin
            for entry in [
                "c-x b[0],b[1],b[2],b[3],q[4]" => [["c-x", "b[0]", "b[1]", "b[2]", "b[3]", "q[4]"]],
                "c-rx b[0],b[1],b[2],q[3],3.14" => [["c-rx", "b[0]", "b[1]", "b[2]", "q[3]", "3.14"]],
                "c-mx90 b[0],b[1],b[2],b[3],q[4]" => [["c-mx90", "b[0]", "b[1]", "b[2]", "b[3]", "q[4]"]],
                "c-sdag b[0],b[1],q[2]" => [["c-sdag", "b[0]", "b[1]", "q[2]"]],
                "c-cnot b[0],b[1],b[2],b[3],b[4],q[5],q[6]" => [["c-cnot", "b[0]", "b[1]", "b[2]", "b[3]", "b[4]", "q[5]", "q[6]"]],
                "c-toffoli b[0],b[1],b[2],b[3],q[4],q[5],q[6]" => [["c-toffoli", "b[0]", "b[1]", "b[2]", "b[3]", "q[4]", "q[5]", "q[6]"]],
                "c-swap b[0],b[1],b[2],q[3],q[4]" => [["c-swap", "b[0]", "b[1]", "b[2]", "q[3]", "q[4]"]],
                "c-crk b[0],b[1],b[2],b[3],b[4],b[5],q[6],q[7],4" => [["c-crk", "b[0]", "b[1]", "b[2]", "b[3]", "b[4]", "b[5]", "q[6]", "q[7]", "4"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end

            for entry in [
                "c-x b[0],bit1,b[2],bit3,q[4]" => [["c-x", "b[0]", "bit1", "b[2]", "bit3", "q[4]"]],
                "c-rx b[0],bit1,b[2],q[3],3.14" => [["c-rx", "b[0]", "bit1", "b[2]", "q[3]", "3.14"]],
                "c-mx90 b[0],bit1,b[2],bit3,q[4]" => [["c-mx90", "b[0]", "bit1", "b[2]", "bit3", "q[4]"]],
                "c-sdag b[0],bit1,q[2]" => [["c-sdag", "b[0]", "bit1", "q[2]"]],
                "c-cnot bit0,b[1],bit2,b[3],b[4],q[5],q[6]" => [["c-cnot", "bit0", "b[1]", "bit2", "b[3]", "b[4]", "q[5]", "q[6]"]],
                "c-toffoli bit0,b[1],bit2,bit3,q[4],q[5],q[6]" => [["c-toffoli", "bit0", "b[1]", "bit2", "bit3", "q[4]", "q[5]", "q[6]"]],
                "c-swap bit0,b[1],bit2,q[3],q[4]" => [["c-swap", "bit0", "b[1]", "bit2", "q[3]", "q[4]"]],
                "c-crk bit0,b[1],bit2,bit3,b[4],bit5,q[6],q[7],4" => [["c-crk", "bit0", "b[1]", "bit2", "bit3", "b[4]", "bit5", "q[6]", "q[7]", "4"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "spaces" begin
            for entry in [
                " swap   qubit9 , q[  2 ] " => [["swap", "qubit9", "q[2]"]],
                "  crk q[9],qubit2,3 " => [["crk", "q[9]", "qubit2", "3"]],
                "c-i bit1 ,    qubit2   " => [["c-i", "bit1", "qubit2"]],
                "    c-toffoli   b[ 0] , q[ 9 ],  q[2  ] , q[ 3 ] " => [["c-toffoli", "b[0]", "q[9]", "q[2]", "q[3]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end

            for entry in [
                " version 1 " => [["version", "1"]],
                "  version  2.3  " => [["version", "2.3"]],
                " qubits  4 " => [["qubits", "4"]],
                "  qubits  52  " => [["qubits", "52"]],
                " map  q[ 0],  data" => [["map", "q[0]", "data"]],
                "  map b[ 1 ] ,bit1" => [["map", "b[1]", "bit1"]],
                " prep_z  q[3 ] " => [["prep_z", "q[3]"]],
                "prep_x q[ 12 ]" => [["prep_x", "q[12]"]],
                "  measure_x  q[4 ]  " => [["measure_x", "q[4]"]],
                "   measure_z    q[ 32]" => [["measure_z", "q[32]"]],
                " measure  q[ 5 ]" => [["measure", "q[5]"]],
                " measure_all " => [["measure_all"]],
                "   measure_parity    q[ 6],  y , q[ 7 ]    , z" => [["measure_parity", "q[6]", "y", "q[7]", "z"]],
                "display" => [["display"]],
                " display  b[ 0 ] " => [["display", "b[0]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "comments" begin
            for entry in [
                "swap   qubit9 , q[  2 ]   # This is a comment   " => [["swap", "qubit9", "q[2]"]],
                "crk q[9],qubit2,3# This is another comment" => [["crk", "q[9]", "qubit2", "3"]],
                "c-i bit1 ,    qubit2 # This a third comment" => [["c-i", "bit1", "qubit2"]],
                "c-toffoli   b[ 0] , q[ 9 ],  q[2  ] , q[ 3 ]#Toomanycommentsinthiscode" => [["c-toffoli", "b[0]", "q[9]", "q[2]", "q[3]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end

            for entry in [
                "version 1 # Comment1" => [["version", "1"]],
                "version  2.3# Comment 2" => [["version", "2.3"]],
                "qubits  4      # Comment 3 and counting" => [["qubits", "4"]],
                "qubits  52  # A bit tired reading sooo much" => [["qubits", "52"]],
                "map  q[ 0],  data#This reminds me of Deep Rock Galactic dwarves comments" => [["map", "q[0]", "data"]],
                "map b[ 1 ] ,bit1#IJUSTDON'TKNOWWHATI'MDOING" => [["map", "b[1]", "bit1"]],
                "prep_z  q[3 ] #JUST cHeckInGG rAandOoOm things" => [["prep_z", "q[3]"]],
                "prep_x q[ 12 ]# I guess this is enough" => [["prep_x", "q[12]"]],
                "display # no, it wasn't" => [["display"]],
                raw"display  b[ 0 ]#Checking weird chars now... #!?$%&'()*+.,:;-/_<>=@[]^`{}|~\\" => [["display", "b[0]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "subcircuits" begin
            for entry in [
                ".init" => [[".init"]],
                ".measure(45)" => [[".measure", "45"]],
                ".grover(3)" => [[".grover", "3"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "Multi Statements" begin
            for entry in [
                "{h q[0]|h q[1]|h q[2]|h q[3]|h oracle}" => [["h", "q[0]", "|", "h", "q[1]", "|", "h", "q[2]", "|", "h", "q[3]", "|", "h", "oracle"]],
                "{ prep_z q[0] | prep_z q[1] | prep_z q[2] }" => [["prep_z", "q[0]", "|", "prep_z", "q[1]", "|", "prep_z", "q[2]"]],
                "{ measure q[0] | measure q[1] | measure q[2] }" => [["measure", "q[0]", "|", "measure", "q[1]", "|", "measure", "q[2]"]],
                "{ h q[0] | prep_z q[1] | measure q[2] }" => [["h", "q[0]", "|", "prep_z", "q[1]", "|", "measure", "q[2]"]],
                "{ x q[0] | x q[1] | x q[2] | x q[3] }" => [["x", "q[0]", "|", "x", "q[1]", "|", "x", "q[2]", "|", "x", "q[3]"]],
                "{ c-i b1,q2 | c-swap b[0],b[10],b[11],q[1],q[2] | toffoli q3,q[6],q[10] }" => [["c-i", "b1", "q2", "|", "c-swap", "b[0]", "b[10]", "b[11]", "q[1]", "q[2]", "|", "toffoli", "q3", "q[6]", "q[10]"]],
                "{  c-i  b1, q2 |  c-swap b[0 ],  b[ 10] ,b[11], q[ 1] ,q[ 2 ]   | toffoli q3,  q[   6 ], q[ 10    ]   }" => [["c-i", "b1", "q2", "|", "c-swap", "b[0]", "b[10]", "b[11]", "q[1]", "q[2]", "|", "toffoli", "q3", "q[6]", "q[10]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end
        end

        @testset "Multiple qubits specification" begin
            for entry in [
                "h q[0:2]" => [["h", "q[0:2]"]],
                "h q[0:2,5]" => [["h", "q[0:2&5]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end

            for entry in [
                "c-i b[0:2],qbit3" => [["c-i", "b[0:2]", "qbit3"]],
                "c-i b[1,4:6,9],q[11]" => [["c-i", "b[1&4:6&9]", "q[11]"]],
                "c-i b[0:2,4:6,9:13],q[20]" => [["c-i", "b[0:2&4:6&9:13]", "q[20]"]],

                "c-i   b[ 0:2 , 4, 7:9], rbit1,  b[3], qbit2" =>[["c-i", "b[0:2&4&7:9]", "rbit1", "b[3]", "qbit2"]],
                "c-swap b[0:2, 6, 8:10], qbit1 , q[2]" => [["c-swap", "b[0:2&6&8:10]", "qbit1", "q[2]"]],
            ]
                @test parseCQASMCode(entry.first) == entry.second
            end

            #  Alternatives for this section of the parser:
            # "h q[0:2]" => [["h", "q[0]", "|", "h", "q[1]", "|", "h", "q[2]"]],
            # "c-i b[0:2],rbit5,b[3],qbit7" => [["c-i", "b[0]", "b[1]", "b[2]", "rbit5", "b[3]", "qbit7"]],
            # "c-i b[0:2,4,6:10,3,14],rbit20,b[17],qbit22" => [["c-i", "b[0]", "b[1]", "b[2]", "b[4]", "b[6]", "b[7]", "b[8]", "b[9]", "b[10]", "b[3]", "b[14]", "rbit20", "b[17]", "qbit22"]],
        end
    end
    # Add complete examples, 1 or 2 will be enough
end

# subcircuits whole examples