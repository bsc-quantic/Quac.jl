module cQASM

export parseCQASMCode

# First imports
import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp;

machine = let
    # Define generic methods
    function addActionToRegEx(strRegEx, enterAction, exitAction)
        strRegEx.actions[:enter] = [enterAction]
        strRegEx.actions[:exit] = [exitAction]
    
        return strRegEx
    end

    function createGates(gate::Symbol, gateName, params, c_gateName, c_params)
        @eval $gate = re.cat($gateName, $params)
        @eval $(Symbol("c_" * String(gate))) = re.cat($c_gateName, $c_params)
    end

    ## Define cQASM syntax ##

    # 1.1 Base strings
    newLine = re"\r?\n"
    zeroOrMoreSpaces = re"[\t ]*"
    oneOrMoreSpaces = re"[\t ]+"
    zeroOrMoreDigits = re"[0-9]*"
    oneOrMoreDigits = re"[0-9]+"
    
    natural = oneOrMoreDigits
    float_custom = re.cat(natural, ".", natural)
    angle = float_custom

    tag = re"[a-z_][a-z0-9_]*"
    idStreamSingle = re.cat(zeroOrMoreSpaces, oneOrMoreDigits, re.opt(re.cat(":", oneOrMoreDigits)), zeroOrMoreSpaces)
    idStreamMultiple = re.cat(idStreamSingle, re.rep(re.cat(",", idStreamSingle)))
    qubitIDStream = re.cat("q[", idStreamMultiple, "]")
    bitIDStream = re.cat("b[", idStreamMultiple, "]")
    qid = re.alt(tag, qubitIDStream)
    bid = re.alt(tag, bitIDStream)

    controlBits = re.rep1(re.cat(bid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces))

    # 1.2 Gates
    # One qubit gates
    for name in [:i, :h, :x, :y, :z, :x90, :y90, :mx90, :my90, :s, :sdag, :t]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :instrEnter, :instrExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid), :paramsEnter, :paramsExit)
        
        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :instrEnter, :instrExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        createGates(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Two qubit gates
    for name in [:cnot, :cz, :swap]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :instrEnter, :instrExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :instrEnter, :instrExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        createGates(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Three qubit gates
    for name in [:toffoli]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :instrEnter, :instrExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :instrEnter, :instrExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        createGates(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Rotation gates
    for name in [:rx, :ry, :rz, :cr]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :instrEnter, :instrExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :instrEnter, :instrExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, zeroOrMoreSpaces,",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle), :paramsEnter, :paramsExit)

        createGates(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Another unique gates
    gateNameRegEx = addActionToRegEx(re.parse("crk"), :instrEnter, :instrExit)
    paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural), :paramsEnter, :paramsExit)

    c_gateNameRegEx = addActionToRegEx(re.parse("c-crk"), :instrEnter, :instrExit)
    c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, bid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural), :paramsEnter, :paramsExit)

    createGates(:crk, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)

    # 1.3 Composed statements
    # version = re.parse("version" * oneOrMoreSpaces * "(" * natural * "|" * float_custom * ")")

    version = re.cat("version", oneOrMoreSpaces, re.alt(natural, float_custom))
    numberQubits = re.cat("qubits", oneOrMoreSpaces, natural)
    #comment = re.parse("#[a-zA-Z0-9_\t ]*")
    comment = re.parse("#[a-zA-Z0-9#[!\"#\$%&'()*+,\\-./:;<=>?@\\[\\]\\^_`{|}~]*\t ]*")
    gate = re.alt(i, h, x, y, z, x90, y90, mx90, my90, s, sdag, t,
                    c_i, c_h, c_x, c_y, c_z, c_x90, c_y90, c_mx90, c_my90, c_s, c_sdag, c_t)
    # gate = re.alt(  i, h, x, y, z, x90, y90, mx90, my90, s, sdag, t,
    #                 c_i, c_h, c_x, c_y, c_z, c_x90, c_y90, c_mx90, c_my90, c_s, c_sdag, c_t,
    #                 cnot, cz, swap,
    #                 c_cnot, c_cz, c_swap,
    #                 toffoli,
    #                 c_toffoli,
    #                 rx, ry, rz, cr,
    #                 c_rx, c_ry, c_rz, c_cr
    #                 crk,
    #                 c_crk)
    # prepState = re"prep_[xyz]" * oneOrMoreSpaces * qid
    # map = re"map" * oneOrMoreSpaces * (bid | qid) * zeroOrMoreSpaces * "," * zeroOrMoreSpaces * tag
    # measurement = re"measure_(?:[xyz]|all)|measure$" * oneOrMoreSpaces * qid
    # display = re"display" * oneOrMoreSpaces * bitName

    # statement = initLine * (version | numberQubits | comment | gate | prepState | map | measurement | display) * endLine
    line = re.cat(re.alt(version, numberQubits, comment, gate), zeroOrMoreSpaces, newLine)
    # line = re.alt(version, numberQubits, comment, oneQubitGate, twoQubitGate, threeQubitGate, rotationGate, uniqueGate) * re.parse(newLine)
    cqasmCode = re.rep(re.alt(line, re.cat(zeroOrMoreSpaces, newLine)))

    # Compile the regex to a FSM
    Automa.compile(cqasmCode)
end;

actions = Dict(
    :instrEnter => :(mark = p),
    :instrExit => :(statementInstr = String(data[mark:p - 1])),
    :paramsEnter => :(mark = p),
    :paramsExit => quote
        statementParams = String(data[mark:p - 1])
        statementParamsSplit = String.(split(replace(statementParams, r"[\t ]" => ""), ","))

        statement = pushfirst!(statementParamsSplit, statementInstr)
        push!(statementsSet, statement)
    end,
);

context = Automa.CodeGenContext();

@eval function parseCQASMCode(data::String)
    data = data * "\n"

    mark = 1
    statementInstr = ""
    statementParams = ""
    statementsSet = Vector{String}[]

    $(Automa.generate_init_code(context, machine))

    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, machine, actions))

    iszero(cs) || error("failed to parse on byte ", p, "\n[last chunk '", data[mark:p - 1], "' from string:\n'", data, "']")
    return statementsSet
end;

parseCQASMCode(String(readchomp("src/Serializers/test1.cq")))


function paint_machine(regexp)
    machine = Automa.compile(regexp)
    write("C:/Users/German/Documents/Work/Quac.jl/src/Serializers/csv.dot", Automa.machine2dot(machine))
end;

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


end