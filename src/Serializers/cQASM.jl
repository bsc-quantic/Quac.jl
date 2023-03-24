module cQASM

export getInstructions

# First imports
import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp;

machine = let
    ## Define cQASM syntax ##

    # 1.1 Base strings
    initLine = "^[\t ]*"
    endLine = "[\t ]*\$"
    newLine = "\r?\n"
    zeroOrMoreSpaces = "[\t ]*"
    oneOrMoreSpaces = "[\t ]+"
    zeroOrMoreDigits = "[0-9]*"
    oneOrMoreDigits = "[0-9]+"
    
    natural = oneOrMoreDigits
    float_custom = natural * "\\." * natural
    angle = float_custom

    tag = "[a-z_][a-z0-9_]*"
    idStreamSingle = zeroOrMoreSpaces * oneOrMoreDigits * "(:" * oneOrMoreDigits * ")?" * zeroOrMoreSpaces
    idStreamMultiple = "(" * idStreamSingle * ")(," * idStreamSingle * ")*"
    qubitIDStream = "q\\[" * idStreamMultiple * "\\]"
    bitIDStream = "b\\[" * idStreamMultiple * "\\]"
    qid = "(" * tag * "|" * qubitIDStream * ")"
    bid = "(" * tag * "|" * bitIDStream * ")"

    controlBits = "(" * bid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * ")+"

    # 1.2 Gates
    # One qubit gates
    for name in [:i, :h, :x, :y, :z, :x90, :y90, :mx90, :my90, :s, :sdag, :t]
        regexStr = re.parse(String(name) * oneOrMoreSpaces * qid)
        @eval $name = $regexStr

        c_regexStr = re.parse("c-" * String(name) * oneOrMoreSpaces * controlBits * zeroOrMoreSpaces * qid)
        @eval $(Symbol("c_" * String(name))) = $c_regexStr
    end

    # Two qubit gates
    for name in [:cnot, :cz, :swap]
        regexStr = re.parse(String(name) * oneOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid)
        @eval $name = $regexStr

        c_regexStr = re.parse("c-" * String(name) * oneOrMoreSpaces * controlBits * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid)
        @eval $(Symbol("c_" * String(name))) = $c_regexStr
    end

    # Three qubit gates
    for name in [:toffoli]
        regexStr = re.parse(String(name) * oneOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid)
        @eval $name = $regexStr

        c_regexStr = re.parse("c-" * String(name) * oneOrMoreSpaces * controlBits * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid)
        @eval $(Symbol("c_" * String(name))) = $c_regexStr
    end

    # Rotation gates
    for name in [:rx, :ry, :rz, :cr]
        regexStr = re.parse(String(name) * oneOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * angle)
        @eval $name = $regexStr

        c_regexStr = re.parse("c-" * String(name) * oneOrMoreSpaces * controlBits * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * angle)
        @eval $(Symbol("c_" * String(name))) = $c_regexStr
    end

    # Another unique gates
    crk = re.parse("crk" * oneOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * natural)
    c_crk = re.parse("c-crk" * oneOrMoreSpaces * bid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * qid * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * natural)

    # 1.3 Composed statements
    version = re.parse("version" * oneOrMoreSpaces * "(" * natural * "|" * float_custom * ")")
    numberQubits = re.parse("qubits" * oneOrMoreSpaces * natural)
    #comment = re.parse("#[a-zA-Z0-9_\t ]*")
    comment = re.parse("#[a-zA-Z0-9#[!\"#\$%&'()*+,\\-./:;<=>?@\\[\\]\\^_`{|}~]*\t ]*")
    oneQubitGate = re.alt(i, h, x, y, z, x90, y90, mx90, my90, s, sdag, t)
    # gate =  i   | h     | x     | y     | z     | rx    | ry    | rz    | x90   | y90   | mx90  | my90  | s     | sdag  | t     | cnot  | toffoli   | cz    | swap  | crk   | cr |
    #         c_i | c_h   | c_x   | c_y   | c_z   | c_rx  | c_ry  | c_rz  | c_x90 | c_y90 | c_mx90| c_my90| c_s   | c_sdag| c_t   | c_cnot| c_toffoli | c_cz  | c_swap| c_crk | c_cr
    # prepState = re"prep_[xyz]" * oneOrMoreSpaces * qid
    # map = re"map" * oneOrMoreSpaces * (bid | qid) * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * tag
    # measurement = re"measure_(?:[xyz]|all)|measure$" * oneOrMoreSpaces * qid
    # display = re"display" * oneOrMoreSpaces * bitName

    # statement = initLine * (version | numberQubits | comment | gate | prepState | map | measurement | display) * endLine
    line = re.alt(version, numberQubits, comment, oneQubitGate) * re.parse(newLine)
    # line = re.alt(version, numberQubits, comment, oneQubitGate, twoQubitGate, threeQubitGate, rotationGate, uniqueGate) * re.parse(newLine)
    cqasmCode = re.rep(line)

    line.actions[:exit] = [:lineExit]
    oneQubitGate.actions[:enter] = [:oneQubitGateEnter]
    oneQubitGate.actions[:exit] = [:oneQubitGateExit]
    # oneQubitGate_qubit.actions[:enter] = [oneQubitGate_qubitEnter]
    # oneQubitGate_qubit.actions[:exit] = [oneQubitGate_qubitExit]

    function paint_machine(regexp)
        machine = Automa.compile(regexp)
        write("C:/Users/German/Documents/Work/Quac.jl/src/Serializers/csv.dot", Automa.machine2dot(machine))
    end;

    # Compile the regex to a FSM
    Automa.compile(cqasmCode)
end;

actions = Dict(
    :lineExit => quote
        statement = String(data[mark:p - 1])
        split(statement)
        push!(statementsSet, statement)
        mark = p
    end,
);

# actions = Dict(
#     :enter => :(mark= p),
#     :exit_header => :(header = String(data[mark+1:p-1])),
#     :exit_seqline => quote
#         doff = length(seqbuffer) + 1
#         resize!(seqbuffer, length(seqbuffer) + p - mark)
#         copyto!(seqbuffer, doff, data, mark, p-mark)
#     end,
#     :exit_record => quote
#         sequence = LongSequence{A}(seqbuffer)
#         empty!(seqbuffer)
#         record = FastaRecord(header, sequence)
#         push!(records, record)
#     end,
# );

# context = Automa.CodeGenContext();
# @eval function parse_fasta(::Type{A}, data::Union{String,Vector{UInt8}}) where {A <: Alphabet}
#     mark = 0
#     records = FastaRecord{A}[]
#     header = ""
#     sequence = LongSequence{A}()
#     seqbuffer = UInt8[]

    
#     (Automa.generate_init_code(context, machine))
#     p_end = p_eof = lastindex(data)
#     (Automa.generate_exec_code(context, machine, actions))

#     iszero(cs) || error("failed to parse on byte ", p)
#     return records
# end;

# g = make_grammar([:expr], flatten(rules, String))

# input = "i q[0]"
# p = parse(g, input)

context = Automa.CodeGenContext();

@eval function parseCQASMCode(data::Union{String,Vector{UInt8}})
    mark = 1
    statementsSet = Array{String[]}; # This should let to push! new arrays inside, probably it is ok already.

    $(Automa.generate_init_code(context, machine))

    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, machine, actions))

    println(statementsSet)

    iszero(cs) || error("failed to parse on byte ", p)
    return nothing
end;

# function parseCQASMCode(code::String)
    # instructions = split(lstrip(rstrip(fileStr)), "\n")
    # for inst in instructions
    #     println("Parsing instruction: " , inst)
    #     parseInstruction(inst)
    # end
# end

end