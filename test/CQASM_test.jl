@testset "Serialization" begin
    @testset "cQASM_generic" begin
        
    end

    @testset "cQASM_gates" begin
        using Quac.cQASM

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

            "c-i b[0],q[0]" => "c-i",
            "c-h b[0],q[0]" => "c-h",
            "c-x b[0],q[0]" => "c-x",
            "c-y b[0],q[0]" => "c-y",
            "c-z b[0],q[0]" => "c-z",
            "c-x90 b[0],q[0]" => "c-x90",
            "c-y90 b[0],q[0]" => "c-y90",
            "c-mx90 b[0],q[0]" => "c-mx90",
            "c-my90 b[0],q[0]" => "c-my90",
            "c-s b[0],q[0]" => "c-s",
            "c-sdag b[0],q[0]" => "c-sdag",
            "c-t b[0],q[0]" => "c-t",
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
end