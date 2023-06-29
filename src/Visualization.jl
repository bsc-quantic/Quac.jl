using Cobweb: h

texname(::Type{Op}) where {Op<:Operator} = String(nameof(Op))

texname(::Type{Sd}) = """S<tspan font-size="60%" baseline-shift="super">†</tspan>"""
texname(::Type{Td}) = """T<tspan font-size="60%" baseline-shift="super">†</tspan>"""

texname(::Type{Rx}) = """R<tspan font-size="60%" baseline-shift="sub">X</tspan>"""
texname(::Type{Ry}) = """R<tspan font-size="60%" baseline-shift="sub">Y</tspan>"""
texname(::Type{Rz}) = """R<tspan font-size="60%" baseline-shift="sub">Z</tspan>"""

texname(::Type{Hz}) = """H<tspan font-size="60%" baseline-shift="sub">Z</tspan>"""

texname(::Type{FSim}) = """F<tspan font-size="60%" baseline-shift="sub">S</tspan>"""

function draw end
export draw

const DRAWING_STYLE = h.style(
    """
    .wire {
        stroke: black;
        stroke-width: 2px;
    }

    .lane {}
    .virtual {}

    .block {
        stroke: black;
        fill: transparent;
    }

    text {
        text-anchor: middle;
        dominant-baseline: central;
    }
""",
    type = "text/css",
)

function __svg_vcat_blocks(blocks...)
    container = h.svg(width = maximum(x -> x.width, blocks), height = sum(x -> parse(Int, x.height), blocks))

    for (block, y) in zip(blocks, Iterators.flatten([0, cumsum(Iterators.map(x -> parse(Int, x.height), blocks))]))
        block.y = y
        push!(container, block)
    end

    return container
end

function __svg_hcat_blocks(blocks...)
    container = h.svg(height = maximum(x -> x.height, blocks), width = sum(x -> parse(Int, x.width), blocks))

    for (block, x) in zip(blocks, Iterators.flatten([0, Iterators.map(x -> parse(Int, x.width), blocks)]))
        block.x = x
        push!(container, block)
    end

    return container
end

function draw(circuit::Circuit; kwargs...)
    n = lanes(circuit)

    if isempty(circuit)
        svg = __svg_vcat_blocks([draw(Gate{I}(lane); kwargs...) for lane in 1:n]...)
        push!(svg, DRAWING_STYLE)
        return svg
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

        return ms
    end |> Iterators.flatten |> collect

    svg = mapreduce(__svg_hcat_blocks, _moments) do moment
        (min, max) = extrema(mapreduce(lanes, ∪, moment))
        moment = [map(I, filter(<(min), 1:n))..., moment..., map(I, filter(>(max), 1:n))...]

        mapreduce(x -> draw(x; kwargs...), __svg_vcat_blocks, moment)
    end

    push!(svg, DRAWING_STYLE)

    return svg
end

Base.show(io::IO, ::MIME"image/svg+xml", circuit::Circuit) = print(io, draw(circuit))

draw(gate::Gate{Op,1,P}) where {Op,P} = draw_block(; top = false, bottom = false)

function draw(gate::Gate{Op,N,P}) where {Op,N,P}
    a, b = extrema(lanes(gate))
    n = b - a + 1
    __svg_vcat_blocks(
        draw_block(; top = true, bottom = false),
        fill(draw_multiblock_mid(), (n - 2))...,
        draw_block(; top = false, bottom = true),
    )
end

draw(::Gate{I,1,NamedTuple{(),Tuple{}}}) =
    h.svg(h.line."wire lane"(; x1 = -25, y1 = 0, x2 = 25, y2 = 0), viewBox = "-25 -25 50 50", width = 50, height = 50)

for Op in [X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Hz, FSim]
    @eval draw(::Gate{$Op,1,P}; kwargs...) where {P} = draw_block(texname($Op); kwargs...)
end

function draw(gate::Gate{<:Control}; kwargs...)
    c = control(gate)
    t = target(gate)
    r = range(extrema(lanes(gate))...)

    # TODO Control{Swap}
    @assert length(t) == 1

    __svg_vcat_blocks(
        [
            if lane == first(r)
                draw_copy(:top)
            elseif lane ∈ c
                draw_copy(:mid)
            else
                draw_cross()
            end for lane in r if lane < only(t)
        ]...,
        draw(
            Gate{targettype(operator(gate))}(target(gate)...; parameters(gate)...);
            top = !any(<(only(t)), c),
            bottom = !any(>(only(t)), c),
            kwargs...,
        ),
        [
            if lane == last(r)
                draw_copy(:bottom)
            elseif lane ∈ c
                draw_copy(:mid)
            else
                draw_cross()
            end for lane in r if lane > only(t)
        ]...,
    )
end

function draw_block(label = ""; top::Bool = false, bottom::Bool = false)
    drawing = h.svg(
        h.line."wire lane"(x1 = -25, y1 = 0, x2 = -15, y2 = 0),
        h.line."wire lane"(x1 = 25, y1 = 0, x2 = 15, y2 = 0),
        h.rect."block"(x = -15, y = -15, width = 30, height = 30),
        h.text."label"(label, x = 0, y = 0),
        viewBox = "-25 -25 50 50",
        width = 50,
        height = 50,
    )
    top && push!(drawing, h.line."wire virtual"(x1 = 0, y1 = 25, x2 = 0, y2 = 15))
    bottom && push!(drawing, h.line."wire virtual"(x1 = 0, y1 = -25, x2 = 0, y2 = -15))
    return drawing
end

draw_multiblock_mid() = h.svg(
    h.line."wire lane"(x1 = -25, y1 = 0, x2 = -15, y2 = 0),
    h.line."wire lane"(x1 = 25, y1 = 0, x2 = 15, y2 = 0),
    h.line(x1 = -25, y1 = 0, x2 = 25, y2 = 0), # TODO assign class. fill?
    h.line(x1 = 0, y1 = -25, x2 = 0, y2 = 25), # TODO assign class. fill?
    viewBox = "-25 -25 50 50",
    width = 50,
    height = 50,
)

draw_cross() = h.svg(
    h.line."wire lane"(x1 = -25, y1 = 0, x2 = 25, y2 = 0),
    h.line."wire virtual"(x1 = 0, y1 = -25, x2 = 0, y2 = 25),
    viewBox = "-25 -25 50 50",
    width = 50,
    height = 50,
)

function draw_copy(dir::Symbol)
    (a, b) = if dir === :top
        0, 25
    elseif dir === :bottom
        0, -25
    elseif dir === :mid
        25, -25
    else
        throw(ArgumentError("`dir`=$dir is invalid"))
    end

    h.svg(
        h.line."wire lane"(x1 = -25, y1 = 0, x2 = 25, y2 = 0),
        h.circle."copy"(cx = 0, cy = 0, r = 5),
        h.line."wire virtual"(x1 = 0, y1 = a, x2 = 0, y2 = b),
        viewBox = "-25 -25 50 50",
        width = 50,
        height = 50,
    )
end
