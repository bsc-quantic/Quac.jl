module cQASM

export parseCQASMCode
# export buildTree

# First imports
import Automa
import Automa.RegExp: @re_str
const re = Automa.RegExp;

# oneQubitGate = re"i|h|x|y|z|x90|y90|mx90|my90|s|sdag|t"
# twoQubitGate = re"cnot|cz|swap"
# threeQubitGate = re"toffoli"
# rotationGate = re"rx|ry|rz|cr"
# uniqueGate = re"crk"
# oneQubitGate_bc = re"c-i|c-h|c-x|c-y|c-z|c-x90|c-y90|c-mx90|c-my90|c-s|c-sdag|c-t"
# twoQubitGate_bc = re"c-cnot|c-cz|c-swap"
# threeQubitGate_bc = re"c-toffoli"
# rotationGate_bc = re"c-rx|c-ry|c-rz|c-cr"
# uniqueGate_bc = re"c-crk"

# keyword = re"version|qubits|map|prep|measure|display"
# identifier = re"[A-Za-z_][0-9A-Za-z_]*"
# comment = re"#[^\r\n]*"
# newline = re"\r?\n"

# minijulia = Automa.compile(
#     re","               => :(emit(:comma)),
#     re":"               => :(emit(:colon)),
#     re"\."              => :(emit(:dot)),
#     re"\("              => :(emit(:lparen)),
#     re"\)"              => :(emit(:rparen)),
#     re"\["              => :(emit(:lbracket)),
#     re"]"               => :(emit(:rbracket)),
#     re"{"               => :(emit(:lbrace)),
#     re"}"               => :(emit(:rbrace)),
#     re"\|"              => :(emit(:vbar)),
#     re"[0-9]+\.[0-9]+"  => :(emit(:float)),
#     re"[0-9]+"          => :(emit(:integer)),
#     re"[\t ]+"          => :(emit(:spaces)),

#     oneQubitGate        => :(emit(:oneQubitGate)),
#     twoQubitGate        => :(emit(:twoQubitGate)),
#     threeQubitGate      => :(emit(:threeQubitGate)),
#     rotationGate        => :(emit(:rotationGate)),
#     uniqueGate          => :(emit(:uniqueGate)),
#     oneQubitGate_bc     => :(emit(:oneQubitGate_bc)),
#     twoQubitGate_bc     => :(emit(:twoQubitGate_bc)),
#     threeQubitGate_bc   => :(emit(:threeQubitGate_bc)),
#     rotationGate_bc     => :(emit(:rotationGate_bc)),
#     uniqueGate_bc       => :(emit(:uniqueGate_bc)),

#     keyword             => :(emit(:keyword)),
#     identifier          => :(emit(:identifier)),
#     comment             => :(emit(:comment)),
#     newline             => :(emit(:newline)),
# )

# context = Automa.CodeGenContext()

# @eval function tokenize(data)
#     $(Automa.generate_init_code(context, minijulia))
#     p_end = p_eof = sizeof(data)
#     tokens = Tuple{Symbol,String}[]
#     emit(kind) = push!(tokens, (kind, data[ts:te]))
#     while p ≤ p_eof && cs > 0
#         $(Automa.generate_exec_code(context, minijulia))
#     end
#     if cs < 0
#         error("failed to tokenize: " * data[ts:te])
#     end
#     return tokens
# end

# a = tokenize(readchomp("src/Serializers/test2.cq"))

# struct TreeState
#     inside_parens::Bool
# end

# function buildTree(data)
#     tokens = tokenize(data)

#     function look_for(::Type{T}, ownType) where T
#         println(tokens[1])
#     end

#     function lookFor(::Type{Float64})

#     end

#     function lookFor(::Type{String})

#     end

#     while !isempty(tokens)
#         token = popfirst!(tokens)
#         println(token)

#         symbol = token[1]
#         string = token[2]

#         if symbol == :comma
    
#         elseif symbol == :colon
#         elseif symbol == :dot
#         elseif symbol == :lparen
#         elseif symbol == :rparen
#         elseif symbol == :lbracket
#         elseif symbol == :rbracket
#         elseif symbol == :lbrace
#             insideBraces = true
#         elseif symbol == :rbrace
#             insideBraces = false
#         elseif symbol == :vbar
#             if insideBraces
#                 # return [reset all variables]
#             end
#         elseif symbol == :integer
#         elseif symbol == :spaces
#             popfirst!(tokens)
#         elseif symbol == :keyword
#             if string == "version "
#                 lookfor(Float64)
#             elseif string == "qubit "
#             end
#         elseif symbol == :identifier
#         elseif symbol == :comment
#         elseif symbol == :newline
#             # return [reset all variables]
#         end

#     end
# end








machine = let
    # Define generic methods
    function addActionToRegEx(strRegEx, enterAction, exitAction=:nothingAction)
        strRegEx.actions[:enter] = [enterAction]
        strRegEx.actions[:exit] = [exitAction]
    
        return strRegEx
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
    
    varName = re"[a-z_][a-z0-9_]*"

    lbracket = addActionToRegEx(re"\[", :lbracketEnter, :lbracketExit)
    rbracket = addActionToRegEx(re"\]", :rbracketEnter, :rbracketExit)
    colon = addActionToRegEx(re":", :colonEnter)
    comma = addActionToRegEx(re",", :commaEnter)
    
    tag = varName
    idStreamSingle = re.cat(zeroOrMoreSpaces, oneOrMoreDigits, re.opt(re.cat(colon, oneOrMoreDigits)), zeroOrMoreSpaces)
    idStreamMultiple = re.cat(idStreamSingle, re.rep(re.cat(comma, idStreamSingle)))
    qubitIDStream = re.cat("q", lbracket, idStreamMultiple, rbracket)
    bitIDStream = re.cat("b", lbracket, idStreamMultiple, rbracket)
    qid = re.alt(tag, qubitIDStream)
    bid = re.alt(tag, bitIDStream)

    lbracket = re"\["
    addActionToRegEx(lbracket, :hey, :heyExit)

    controlBits = re.rep1(re.cat(bid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces))

    # 1.2 Gates
    # One qubit gates
    one_qubit_gates_names = ["i", "h", "x", "y", "z", "x90", "y90", "mx90", "my90", "s", "sdag", "t"]
    
    names_non_controlled_re = re.parse(join(one_qubit_gates_names, "|"))
    names_bit_controlled_re = re.parse(join(["c-" * gate for gate in one_qubit_gates_names], "|"))
    params_non_controlled_re = re.cat(oneOrMoreSpaces, qid)
    params_bit_controlled_re = re.cat(oneOrMoreSpaces, controlBits, qid)

    one_qubit_gates_names_non_controlled = addActionToRegEx(names_non_controlled_re, :statementEnter, :statementExit)
    one_qubit_gates_names_bit_controlled = addActionToRegEx(names_bit_controlled_re, :statementEnter, :statementExit)
    one_qubit_gates_params_non_controlled = addActionToRegEx(params_non_controlled_re, :paramsEnter, :paramsExit)
    one_qubit_gates_params_bit_controlled = addActionToRegEx(params_bit_controlled_re, :paramsEnter, :paramsExit)

    one_qubit_gates_regExp_non_controlled = re.cat(one_qubit_gates_names_non_controlled, one_qubit_gates_params_non_controlled)
    one_qubit_gates_regExp_bit_controlled = re.cat(one_qubit_gates_names_bit_controlled, one_qubit_gates_params_bit_controlled)

    # Two qubit gates
    two_qubit_gates_names = ["cnot", "cz", "swap"]

    names_non_controlled_re = re.parse(join(two_qubit_gates_names, "|"))
    names_bit_controlled_re = re.parse(join(["c-" * gate for gate in two_qubit_gates_names], "|"))
    params_non_controlled_re = re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid)
    params_bit_controlled_re = re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid)

    two_qubit_gates_names_non_controlled = addActionToRegEx(names_non_controlled_re, :statementEnter, :statementExit)
    two_qubit_gates_names_bit_controlled = addActionToRegEx(names_bit_controlled_re, :statementEnter, :statementExit)
    two_qubit_gates_params_non_controlled = addActionToRegEx(params_non_controlled_re, :paramsEnter, :paramsExit)
    two_qubit_gates_params_bit_controlled = addActionToRegEx(params_bit_controlled_re, :paramsEnter, :paramsExit)

    two_qubit_gates_regExp_non_controlled = re.cat(two_qubit_gates_names_non_controlled, two_qubit_gates_params_non_controlled)
    two_qubit_gates_regExp_bit_controlled = re.cat(two_qubit_gates_names_bit_controlled, two_qubit_gates_params_bit_controlled)

    # Three qubit gates
    three_qubit_gates_names = ["toffoli"]

    names_non_controlled_re = re.parse(join(three_qubit_gates_names, "|"))
    names_bit_controlled_re = re.parse(join(["c-" * gate for gate in three_qubit_gates_names], "|"))
    params_non_controlled_re = re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid)
    params_bit_controlled_re = re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid)

    three_qubit_gates_names_non_controlled = addActionToRegEx(names_non_controlled_re, :statementEnter, :statementExit)
    three_qubit_gates_names_bit_controlled = addActionToRegEx(names_bit_controlled_re, :statementEnter, :statementExit)
    three_qubit_gates_params_non_controlled = addActionToRegEx(params_non_controlled_re, :paramsEnter, :paramsExit)
    three_qubit_gates_params_bit_controlled = addActionToRegEx(params_bit_controlled_re, :paramsEnter, :paramsExit)

    three_qubit_gates_regExp_non_controlled = re.cat(three_qubit_gates_names_non_controlled, three_qubit_gates_params_non_controlled)
    three_qubit_gates_regExp_bit_controlled = re.cat(three_qubit_gates_names_bit_controlled, three_qubit_gates_params_bit_controlled)

    # Rotation gates
    rotation_gates_names = ["rx", "ry", "rz", "cr"]

    names_non_controlled_re = re.parse(join(rotation_gates_names, "|"))
    names_bit_controlled_re = re.parse(join(["c-" * gate for gate in rotation_gates_names], "|"))
    params_non_controlled_re = re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle)
    params_bit_controlled_re = re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle)

    rotation_gates_names_non_controlled = addActionToRegEx(names_non_controlled_re, :statementEnter, :statementExit)
    rotation_gates_names_bit_controlled = addActionToRegEx(names_bit_controlled_re, :statementEnter, :statementExit)
    rotation_gates_params_non_controlled = addActionToRegEx(params_non_controlled_re, :paramsEnter, :paramsExit)
    rotation_gates_params_bit_controlled = addActionToRegEx(params_bit_controlled_re, :paramsEnter, :paramsExit)

    rotation_gates_regExp_non_controlled = re.cat(rotation_gates_names_non_controlled, rotation_gates_params_non_controlled)
    rotation_gates_regExp_bit_controlled = re.cat(rotation_gates_names_bit_controlled, rotation_gates_params_bit_controlled)

    # Another unique gates (crk)
    unique_gates_names = ["crk"]

    names_non_controlled_re = re.parse(join(unique_gates_names, "|"))
    names_bit_controlled_re = re.parse(join(["c-" * gate for gate in unique_gates_names], "|"))
    params_non_controlled_re = re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural)
    params_bit_controlled_re = re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural)

    unique_gates_names_non_controlled = addActionToRegEx(names_non_controlled_re, :statementEnter, :statementExit)
    unique_gates_names_bit_controlled = addActionToRegEx(names_bit_controlled_re, :statementEnter, :statementExit)
    unique_gates_params_non_controlled = addActionToRegEx(params_non_controlled_re, :paramsEnter, :paramsExit)
    unique_gates_params_bit_controlled = addActionToRegEx(params_bit_controlled_re, :paramsEnter, :paramsExit)

    unique_gates_regExp_non_controlled = re.cat(unique_gates_names_non_controlled, unique_gates_params_non_controlled)
    unique_gates_regExp_bit_controlled = re.cat(unique_gates_names_bit_controlled, unique_gates_params_bit_controlled)

    # 1.3 Misc statements
    versionName = addActionToRegEx(re.parse("version"), :statementEnter, :statementExit)
    versionParams = addActionToRegEx(re.cat(oneOrMoreSpaces, re.alt(natural, float_custom)), :paramsEnter, :paramsExit)
    version = re.cat(versionName, versionParams)

    numberQubitsName = addActionToRegEx(re.parse("qubits"), :statementEnter, :statementExit)
    numberQubitsParams = addActionToRegEx(re.cat(oneOrMoreSpaces, natural), :paramsEnter, :paramsExit)
    numberQubits = re.cat(numberQubitsName, numberQubitsParams)
    
    comment = re.parse(raw"#[a-zA-Z0-9#!?$%&'()*+.,:;\-/\\_<>=@\[\]^`{}|~ ]*")    # ToDo as it is now, comments does not include the chars ["¡¿´], and maybe some more...

    mappingName = addActionToRegEx(re.parse("map"), :statementEnter, :statementExit)
    mappingParams = addActionToRegEx(re.cat(oneOrMoreSpaces, re.alt(bid, qid), zeroOrMoreSpaces, ",", zeroOrMoreSpaces, tag), :paramsEnter, :paramsExit)
    mapping = re.cat(mappingName, mappingParams)

    prepName = addActionToRegEx(re.parse("prep_[xyz]"), :statementEnter, :statementExit)
    prepParams = addActionToRegEx(re.cat(oneOrMoreSpaces, qid), :paramsEnter, :paramsExit)
    prep = re.cat(prepName, prepParams)

    measureXYZName = addActionToRegEx(re.cat("measure", re.opt(re"_[xyz]")), :statementEnter, :statementExit)
    measureXYZParams = addActionToRegEx(re.cat(oneOrMoreSpaces, qid), :paramsEnter, :paramsExit)
    measureXYZ = re.cat(measureXYZName, measureXYZParams)

    measureAll = addActionToRegEx(re.parse("measure_all"), :statementEnter, :measureAllExit)

    measure_parityName = addActionToRegEx(re.parse("measure_parity"), :statementEnter, :statementExit)
    measure_parityParams = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, re"[xyz]", zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, re"[xyz]"), :paramsEnter, :paramsExit)
    measure_parity = re.cat(measure_parityName, measure_parityParams)

    measure = re.alt(measureXYZ, measureAll, measure_parity)

    display = addActionToRegEx(re.cat("display", re.opt(re.cat(oneOrMoreSpaces, bid))), :statementEnter, :displayExit)

    subCircuitNumberIters = addActionToRegEx(re"[0-9]+", :numberItersEnter, :numberItersExit)
    subcircuit = addActionToRegEx(re.cat(".", varName, re.opt(re.cat("(", subCircuitNumberIters, ")"))), :statementEnter, :subcircuitExit)

    multistatement = addActionToRegEx(re"{[a-zA-Z0-9\t |\-_\[\],]*}", :multiStatementEnter, :multiStatementExit)

    gate = re.alt(  one_qubit_gates_regExp_non_controlled, one_qubit_gates_regExp_bit_controlled,
                    two_qubit_gates_regExp_non_controlled, two_qubit_gates_regExp_bit_controlled,
                    three_qubit_gates_regExp_non_controlled, three_qubit_gates_regExp_bit_controlled,
                    rotation_gates_regExp_non_controlled, rotation_gates_regExp_bit_controlled,
                    unique_gates_regExp_non_controlled, unique_gates_regExp_bit_controlled)

    allStatements = re.cat(zeroOrMoreSpaces, re.alt(version, numberQubits, comment, mapping, prep, measure, display, subcircuit, multistatement, gate), zeroOrMoreSpaces)

    line = re.cat(allStatements, re.opt(comment), newLine)
    cqasmCode = re.rep(re.alt(line, re.cat(zeroOrMoreSpaces, newLine)))

    # Compile the regex to a FSM
    Automa.compile(cqasmCode)
end;

function fixRanges(statementParams)
    statementParamsFixed = ""

    insideBrackets = false
    for char in statementParams
        if char == '[' insideBrackets = true
        elseif char == ']' insideBrackets = false
        end

        if char == ',' && insideBrackets statementParamsFixed *= '&'
        else statementParamsFixed *= char
        end
    end
    
    return statementParamsFixed
end

actions = Dict(
    :nothingAction => :(),

    :lbracketEnter => :(),
    :lbracketExit => :(),
    :rbracketEnter => :(),
    :rbracketExit => :(),
    :colonEnter => :(),
    :commaEnter => :(),

    :statementEnter => :(mark = p),
    :statementExit => :(statementInstr = String(data[mark:p - 1])),

    :paramsEnter => :(mark = p),
    :paramsExit => quote
        statementParams = String(data[mark:p - 1])

        statementParamsFixed = fixRanges(statementParams)                       # Needed to keep multiple qubits specification (e.g. q[0:2,4,7])
        statementParamsSplit = String.(split(replace(statementParamsFixed, r"[\t ]" => ""), ","))

        statementArray = pushfirst!(statementParamsSplit, statementInstr)
        push!(statementsSet, statementArray)
    end,
    
    :measureAllExit => :(push!(statementsSet, ["measure_all"])),
    :displayExit => quote
        statement = String(data[mark:p - 1])
        statementTrimmed = replace(statement, r"[\t ]" => "")
        if length(statementTrimmed) == length("display")
            push!(statementsSet, ["display"])
        else
            pop!(statementsSet)     # This is bc a FSM cannot look forward, so it cannot handle multiple dispatch (the 'display' command is MD). So the FSM needed to add the "display" keyword as it encountered it along the flow. Now, if that statement have parameters, you should rm the previously added entry and add the actual one with parameters.
            push!(statementsSet, ["display", statementTrimmed[length("display") + 1:end]])
        end
    end,
    :numberItersEnter => :(mark2 = p),
    :numberItersExit => :(subStatement = data[mark2:p - 1]),
    :subcircuitExit => quote
        statementInstr = String(data[mark:p - 1])
        statementInstr = split(statementInstr, "(")[1]

        if length(subStatement) > 0
            push!(statementsSet, [statementInstr, subStatement])
        else
            push!(statementsSet, [statementInstr])
        end
    end,
    :multiStatementEnter => :(mark = p),
    :multiStatementExit => quote
        multiStatement = String(data[mark:p - 1])
        push!(statementsSet, [multiStatement])
    end,
);

context = Automa.CodeGenContext();

@eval function parseSingleStatement(data::String)
    mark = mark2 = 1
    statementInstr = statementParams = subStatement = ""
    statementsSet = Vector{String}[]

    $(Automa.generate_init_code(context, machine))

    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, machine, actions))

    iszero(cs) || error("failed to parse on byte ", p, "  [last chunk '", data[mark:p - 1], "' from string '", data, "']")
    return statementsSet
end

function parseMultiStatement(statementsSet::Vector{Array{String, 1}})
    resultingSet = Vector{String}[]
    for statementVector in statementsSet
        statement = statementVector[1]
        if statement[1] == '{'
            parallelStatements = replace(statement, r"[{}]" => "", "|" => "\n")
            multiStatementsSet = parseSingleStatement(parallelStatements * "\n")

            barSeparatedStatement = reduce(vcat, [vcat(statement, "|") for statement in multiStatementsSet])
            push!(resultingSet, barSeparatedStatement[1:length(barSeparatedStatement) - 1])
        else
            push!(resultingSet, statementVector)
        end
    end

    return resultingSet
end

function parseCQASMCode(data::String)
    data = lowercase(data * "\n")

    statementsSet = parseSingleStatement(data)
    statementsSet = parseMultiStatement(statementsSet)

    return statementsSet
end;



a = parseCQASMCode(String(readchomp("src/Serializers/test1.cq")))

function paint_machine(regexp::Automa.RegExp.RE)
    machine = Automa.compile(regexp)
    write("C:/Users/German/Documents/Work/Quac.jl/src/Serializers/csv1.dot", Automa.machine2dot(machine))
end;

function paint_machine(machine::Automa.Machine)
    write("C:/Users/German/Documents/Work/Quac.jl/src/Serializers/csv1.dot", Automa.machine2dot(machine))
end;

end