using Luxor
using MathTeXEngine
using LaTeXStrings

texname(::Type{Op}) where {Op<:Operator} = LaTeXString(String(nameof(Op)))

texname(::Type{Sd}) = L"S^\dagger"
texname(::Type{Td}) = L"T^\dagger"

texname(::Type{Rx}) = L"R_X"
texname(::Type{Ry}) = L"R_Y"
texname(::Type{Rz}) = L"R_Z"

function draw end
export draw

function draw_circuit_lines(num_lanes::Int, num_moments::Int)
    @drawsvg begin
        origin()
        for lane in 1:num_lanes
            line(Point(-25, 50 * (lane - 1)), Point(50 * num_moments - 25, 50 * (lane - 1)), action = :stroke)
        end
    end 50 50 * num_lanes
end

function draw(circuit::Circuit; kwargs...)
    num_qubits = lanes(circuit)

    if isempty(circuit)
        return vcat([draw(Gate{I}(lane); kwargs...) for lane in 1:num_qubits]...)
    end

    # split moments if gates overlap in 1D
    _moments = Iterators.map(moments(circuit)) do moment
        queue = copy(moment)
        ms = Vector{Gate}[]

        # group gates with disjoint lane ranges
        while !isempty(queue)
            ref = popfirst!(queue)
            refrange = range(extrema(lanes(ref))...)
            moment = Gate[ref]
            for gate in queue
                if isdisjoint(refrange, range(extrema(lanes(gate))...))
                    push!(moment, gate)
                end
            end

            queue = setdiff(queue, moment)
            push!(ms, moment)
        end

        ms
    end |> Iterators.flatten |> collect

    @drawsvg begin
        origin()
        draw_circuit(circuit, num_qubits, _moments; kwargs...)
    end 50 * (2 * length(_moments) + 1) 50 * (num_qubits + 1)
end



function draw_circuit(circuit::Circuit, num_qubits, _moments; background = nothing, kwargs...)
    # Draw background if specified
    if background !== nothing
        sethue(background)
        box(Point(0, 0), 25 * (2 * length(circuit) + 1) + 50, 50 * (num_qubits + 1), :fill)
    end

    # Draw horizontal lines for all qubits
    sethue("black")
    translate(-25, 0) # Adjust the starting position
    for i in 1:num_qubits+1
        line(Point(-25, 50 * (i - 1)), Point(25 * (2 * length(circuit) + 1), 50 * (i - 1)), action=:stroke)
    end

    # Draw gates
    for (col, moment) in enumerate(_moments)
        for gate in moment
            x_offset = 50 * (2 * col - 1)
            for lane in lanes(gate)
                y_offset = 50 * (lane - 1)

                # Call the corresponding shape drawing function for the gate
                if gate isa Gate{I, 1, NamedTuple{(), Tuple{}}}
                    # Do nothing for I gate
                elseif gate isa Gate{<:Control}
                    if lane == first(control(gate))
                        draw_copy_shape(:top, x_offset, y_offset; kwargs...)
                    elseif lane ∈ control(gate)[2:end]
                        draw_copy_shape(:mid, x_offset, y_offset; kwargs...)
                    else
                        draw_cross_shape(x_offset, y_offset; kwargs...)
                    end
                else
                    draw_block_shape(texname(operator(gate)); x_offset, y_offset, kwargs...)
                end
            end
        end
    end
end





# Define the draw_moment function to place gates along the lines
function draw_moment(moment, idx; kwargs...)

    (min, max) = extrema(mapreduce(lanes, ∪, moment))
    moment = [map(I, filter(<(min), 1:n))..., moment..., map(I, filter(>(max), 1:n))...]

    vcat(map((gate, lane) -> draw(gate; x_offset = 50 * (idx - 1), y_offset = 50 * (lane - 1), kwargs...), moment, 1:n)...)
end


function Base.show(io::IO, ::MIME"image/svg+xml", circuit::Circuit)
    _ = draw(circuit)
    print(io, svgstring())
end

function draw(gate::Gate{Op,1,P}; x_offset = 0, y_offset = 0, kwargs...) where {Op,P}
    draw_block(x_offset = x_offset, y_offset = y_offset; top = false, bottom = false, kwargs...)
end

function draw(gate::Gate{Op,N,P}; x_offset = 0, y_offset = 0, kwargs...) where {Op,N,P}
    a, b = extrema(lanes(gate))
    n = b - a + 1
    vcat(
        draw_block(x_offset = x_offset, y_offset = y_offset; top = true, bottom = false, kwargs...),
        fill(draw_multiblock_mid(x_offset, y_offset; kwargs...), (n - 2))...,
        draw_block(x_offset = x_offset, y_offset = y_offset; top = false, bottom = true, kwargs...),
    )
end

function draw(::Gate{I,1,NamedTuple{(),Tuple{}}}; x_offset = 0, y_offset = 0, background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        translate(x_offset, y_offset)
        origin()
    end 50 50
end

for Op in [X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz]
    @eval draw(::Gate{$Op,1,P}; x_offset = 0, y_offset = 0, kwargs...) where {P} = draw_block(texname($Op); x_offset = x_offset, y_offset = y_offset, kwargs...)
end

function draw(gate::Gate{<:Control}; x_offset = 0, y_offset = 0, kwargs...)
    c = control(gate)
    t = target(gate)
    r = range(extrema(lanes(gate))...)

    # TODO Control{Swap}
    @assert length(t) == 1

    vcat(
        [
            if lane == first(r)
                draw_copy(:top, x_offset, y_offset; kwargs...)
            elseif lane ∈ c
                draw_copy(:mid, x_offset, y_offset; kwargs...)
            else
                draw_cross(x_offset, y_offset; kwargs...)
            end for lane in r if lane < only(t)
        ]...,
        draw(
            Gate{targettype(operator(gate))}(target(gate)...; parameters(gate)...);
            x_offset = x_offset, y_offset = y_offset,
            top = !any(<(only(t)), c),
            bottom = !any(>(only(t)), c),
            kwargs...,
        ),
        [
            if lane == last(r)
                draw_copy(:bottom, x_offset, y_offset; kwargs...)
            elseif lane ∈ c
                draw_copy(:mid, x_offset, y_offset; kwargs...)
            else
                draw_cross(x_offset, y_offset; kwargs...)
            end for lane in r if lane > only(t)
        ]...,
    )
end

function draw_block(label = ""; x_offset = 0, y_offset = 0, top::Bool = false, bottom::Bool = false, background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        translate(x_offset, y_offset)
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action = :stroke)
        line(Point(25, 0), Point(15, 0), action = :stroke)

        # control connectors
        if top
            line(Point(0, 25), Point(0, 15), action = :stroke)
        end
        if bottom
            line(Point(0, -25), Point(0, -15), action = :stroke)
        end

        rect(-15, -15, 30, 30, action = :stroke)

        # label
        fontsize(16)
        text(label, Point(0, 0), valign = :middle, halign = :center)
    end 50 50
end

function draw_multiblock_mid(x_offset = 0, y_offset = 0; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        translate(x_offset, y_offset)
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action = :stroke)
        line(Point(25, 0), Point(15, 0), action = :stroke)

        # vertical lines
        line(Point(-15, -25), Point(-15, 25), action = :stroke)
        line(Point(15, -25), Point(15, 25), action = :stroke)
    end 50 50
end

function draw_cross(x_offset = 0, y_offset = 0; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        translate(x_offset, y_offset)
        origin()

        line(Point(-25, 0), Point(25, 0), action = :stroke)
        line(Point(0, -25), Point(0, 25), action = :stroke)
    end 50 50
end
function draw_block_shape(label = ""; x_offset = 0, y_offset = 0, top::Bool = false, bottom::Bool = false, kwargs...)
    # lane wire
    line(Point(-25 + x_offset, y_offset), Point(-15 + x_offset, y_offset), action = :stroke)
    line(Point(25 + x_offset, y_offset), Point(15 + x_offset, y_offset), action = :stroke)

    # control connectors
    if top
        line(Point(x_offset, 25 + y_offset), Point(x_offset, 15 + y_offset), action = :stroke)
    end
    if bottom
        line(Point(x_offset, -25 + y_offset), Point(x_offset, -15 + y_offset), action = :stroke)
    end

    rect(-15 + x_offset, -15 + y_offset, 30, 30, action = :stroke)

    # label
    fontsize(16)
    text(label, Point(x_offset, y_offset), valign = :middle, halign = :center)
end

function draw_multiblock_mid_shape(x_offset = 0, y_offset = 0; kwargs...)
    # lane wire
    line(Point(-25 + x_offset, y_offset), Point(-15 + x_offset, y_offset), action = :stroke)
    line(Point(25 + x_offset, y_offset), Point(15 + x_offset, y_offset), action = :stroke)

    # vertical lines
    line(Point(-15 + x_offset, -25 + y_offset), Point(-15 + x_offset, 25 + y_offset), action = :stroke)
    line(Point(15 + x_offset, -25 + y_offset), Point(15 + x_offset, 25 + y_offset), action = :stroke)
end

function draw_cross_shape(x_offset = 0, y_offset = 0; kwargs...)
    line(Point(-25 + x_offset, y_offset), Point(25 + x_offset, y_offset), action = :stroke)
    line(Point(x_offset, -25 + y_offset), Point(x_offset, 25 + y_offset), action = :stroke)
end

function draw_copy_shape(dir::Symbol, x_offset = 0, y_offset = 0; kwargs...)
    line(Point(-25 + x_offset, y_offset), Point(25 + x_offset, y_offset), action = :stroke)
    circle(x_offset, y_offset, 5, action = :fill)

    (a, b) = if dir == :top
        Point(x_offset, y_offset), Point(x_offset, 25 + y_offset)
    elseif dir == :bottom
        Point(x_offset, y_offset), Point(x_offset, -25 + y_offset)
    elseif dir == :mid
        Point(x_offset, 25 + y_offset), Point(x_offset, -25 + y_offset)
    else
        throw(ArgumentError("`dir`=$dir is invalid"))
    end
    line(a, b, action = :stroke)
end


function draw_copy(dir::Symbol, x_offset = 0, y_offset = 0; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        translate(x_offset, y_offset)
        origin()
        line(Point(-25, 0), Point(25, 0), action = :stroke)

        circle(0, 0, 5, action = :fill)

        (a, b) = if dir == :top
            Point(0, 0), Point(0, 25)
        elseif dir == :bottom
            Point(0, 0), Point(0, -25)
        elseif dir == :mid
            Point(0, 25), Point(0, -25)
        else
            throw(ArgumentError("`dir`=$dir is invalid"))
        end
        line(a, b, action = :stroke)
    end 50 50
end