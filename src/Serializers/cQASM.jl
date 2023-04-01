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

    function createStatement(statement::Symbol, statementName, params)
        @eval $statement = re.cat($statementName, $params)
    end

    function createStatement(statement::Symbol, statementName, params, c_statementName, c_params)
        @eval $statement = re.cat($statementName, $params)
        @eval $(Symbol("c_" * String(statement))) = re.cat($c_statementName, $c_params)
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
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :statementEnter, :statementExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid), :paramsEnter, :paramsExit)
        
        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :statementEnter, :statementExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, qid), :paramsEnter, :paramsExit)

        createStatement(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Two qubit gates
    for name in [:cnot, :cz, :swap]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :statementEnter, :statementExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :statementEnter, :statementExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        createStatement(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Three qubit gates
    for name in [:toffoli]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :statementEnter, :statementExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :statementEnter, :statementExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid), :paramsEnter, :paramsExit)

        createStatement(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Rotation gates
    for name in [:rx, :ry, :rz, :cr]
        gateNameRegEx = addActionToRegEx(re.parse(String(name)), :statementEnter, :statementExit)
        paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle), :paramsEnter, :paramsExit)

        c_gateNameRegEx = addActionToRegEx(re.parse("c-" * String(name)), :statementEnter, :statementExit)
        c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, angle), :paramsEnter, :paramsExit)

        createStatement(name, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)
    end

    # Another unique gates
    gateNameRegEx = addActionToRegEx(re.parse("crk"), :statementEnter, :statementExit)
    paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural), :paramsEnter, :paramsExit)

    c_gateNameRegEx = addActionToRegEx(re.parse("c-crk"), :statementEnter, :statementExit)
    c_paramsRegEx = addActionToRegEx(re.cat(oneOrMoreSpaces, controlBits, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, qid, zeroOrMoreSpaces, ",", zeroOrMoreSpaces, natural), :paramsEnter, :paramsExit)

    createStatement(:crk, gateNameRegEx, paramsRegEx, c_gateNameRegEx, c_paramsRegEx)

    # 1.3 Misc statements
    versionName = addActionToRegEx(re.parse("version"), :statementEnter, :statementExit)
    versionParams = addActionToRegEx(re.cat(oneOrMoreSpaces, re.alt(natural, float_custom)), :paramsEnter, :paramsExit)
    version = re.cat(versionName, versionParams)

    numberQubitsName = addActionToRegEx(re.parse("qubits"), :statementEnter, :statementExit)
    numberQubitsParams = addActionToRegEx(re.cat(oneOrMoreSpaces, natural), :paramsEnter, :paramsExit)
    numberQubits = re.cat(numberQubitsName, numberQubitsParams)
    
    comment = re.parse(raw"#[a-zA-Z0-9#!¡?$%&'()*+.,:;\-/\\_<>=@\[\]^`´{}|~ ]*")    # ToDo: as it is now, comments does not include the chars " and ¿, and maybe some more...

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

    # displayAll = addActionToRegEx(re.parse("display"), :statementEnter, :displayAllExit)

    display = addActionToRegEx(re.cat("display", re.opt(re.cat(oneOrMoreSpaces, bid))), :statementEnter, :displayExit)

    # displayBitName = addActionToRegEx(re.parse("display"), :statementEnter, :statementExit)
    # displayBitParams = addActionToRegEx(re.cat(oneOrMoreSpaces, bid), :paramsEnter, :paramsExit)
    # displayBit = re.cat(displayBitName, displayBitParams)

    # display = re.alt(displayAll, displayBit)

    gate = re.alt(i, h, x, y, z, x90, y90, mx90, my90, s, sdag, t,
                    c_i, c_h, c_x, c_y, c_z, c_x90, c_y90, c_mx90, c_my90, c_s, c_sdag, c_t,
                    cnot, cz, swap,
                    c_cnot, c_cz, c_swap,
                    toffoli,
                    c_toffoli,
                    rx, ry, rz, cr,
                    c_rx, c_ry, c_rz, c_cr,
                    crk,
                    c_crk)

    line = re.cat(re.alt(display), zeroOrMoreSpaces, re.opt(comment), newLine)
    cqasmCode = re.rep(re.alt(line, re.cat(zeroOrMoreSpaces, newLine)))

    # Compile the regex to a FSM
    Automa.compile(cqasmCode)
end;

actions = Dict(
    :statementEnter => :(mark = p),
    :statementExit => :(statementInstr = String(data[mark:p - 1])),

    :paramsEnter => :(mark = p),
    :paramsExit => quote
        statementParams = String(data[mark:p - 1])
        statementParamsSplit = String.(split(replace(statementParams, r"[\t ]" => ""), ","))

        statementArray = pushfirst!(statementParamsSplit, statementInstr)
        push!(statementsSet, statementArray)
    end,
    
    :measureAllExit => :(push!(statementsSet, ["measure_all"])),
    :displayExit => quote
        statement = String(data[mark:p - 1])
        statementSplit = String.(split(statement))

        if length(statementSplit) > 1 pop!(statementsSet) end   # This is bc a FSM cannot look forward, so it cannot handle multiple dispatch (the 'display' command is MD). So the FSM needed to add the "display" keyword as it encountered it along the flow. Now, if that statement have parameters, you should rm the previously added entry and add the actual one with parameters.

        push!(statementsSet, statementSplit)
    end,
);

context = Automa.CodeGenContext();

@eval function parseCQASMCode(data::String)
    data = lowercase(data * "\n")

    mark = 1
    statementInstr = ""
    statementParams = ""
    statementsSet = Vector{String}[]

    $(Automa.generate_init_code(context, machine))

    p_end = p_eof = lastindex(data)

    $(Automa.generate_exec_code(context, machine, actions))

    iszero(cs) || return string("failed to parse on byte ", p, "  [last chunk '", data[mark:p - 1], "' from string: '", data, "']")
    return statementsSet
end;

parseCQASMCode(String(readchomp("src/Serializers/test2.cq")))


function paint_machine(regexp)
    machine = Automa.compile(regexp)
    write("C:/Users/German/Documents/Work/Quac.jl/src/Serializers/csv.dot", Automa.machine2dot(machine))
end;

end