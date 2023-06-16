function parse_qflex(filename; sites = nothing)
    open(filename, "r") do io
        n = parse(Int, strip(readline(io)))
        n > 0 || throw(ErrorException("number of qubits must be positive"))

        circuit = Circuit(n)

        # qflex qubit id to Quac qubit number
        mapping = Dict(splat(Pair).(reverse.(enumerate(sort(collect(isnothing(sites) ? (1:n) : sites))))))

        for line in readlines(io)
            # remove trailing spaces, newline characters and moment number
            line = lstrip(isnumeric, strip(line))
            isempty(line) && continue

            gate = if (capture = match(r"x_1_2 (?<target>\d+)", line); !isnothing(capture))
                Rx(mapping[capture["target"]]; θ = π / 2)
            elseif (capture = match(r"y_1_2 (?<target>\d+)", line); !isnothing(capture))
                Ry(mapping[capture["target"]]; θ = π / 2)
            elseif (capture = match(r"hz_1_2 (?<target>\d+)", line); !isnothing(capture))
                # NOTE assume constant parameters for Hz gate
                Hz(mapping[capture["target"]]; θ = π / 4, ϕ = π / 2) # TODO ?
            elseif (capture = match(r"rz\((?<θ>[-]\d+\.\d+)\) (?<target>\d+)", line); !isnothing(capture))
                Rz(mapping[capture["target"]]; θ = capture["θ"])
            elseif (
                capture =
                    match(r"fsim\((?<θ>[-]\d+\.\d+),[ ](?<ϕ>[-]\d+\.\d+)\) (?<source>\d+) (?<target>\d+)", line);
                !isnothing(capture)
            )
                FSim(mapping[capture["source"]], mapping[capture["target"]]; θ = capture["θ"], ϕ = capture["ϕ"])
            else
                error(ErrorException("invalid code: $line"))
            end
            push!(circuit, gate)
        end
    end

    return circuit
end
