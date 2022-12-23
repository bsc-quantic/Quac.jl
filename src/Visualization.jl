using Luxor
using MathTeXEngine

function draw end
export draw

function draw(circuit::Circuit)
    n = lanes(circuit)

    if isempty(circuit)
        return vcat([draw(I) for _ in 1:n]...)
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

    return mapreduce(hcat, _moments) do moment
        (min, max) = extrema(mapreduce(lanes, ∪, moment))
        moment = [map(I, filter(<(min), 1:n))..., moment..., map(I, filter(>(max), 1:n))...]

        mapreduce(draw, vcat, moment)
    end
end

function Base.show(io::IO, ::MIME"image/svg+xml", circuit::Circuit)
    _ = draw(circuit)
    print(io, svgstring())
end

function draw(gate::Gate; top::Bool = false, bottom::Bool = false)
    n = (length ∘ lanes)(gate)
    if n == 1
        draw_block(; top = top, bottom = bottom)
    else
        a, b = extrema(lanes(gate))
        n = b - a + 1
        vcat(draw_block(; top = top), fill(draw_multiblock_mid(), (n - 2))..., draw_block(; bottom = bottom))
    end
end

draw(::I) = draw(I)
function draw(::Type{I})
    @drawsvg begin
        origin()
        line(Point(-25, 0), Point(25, 0), action = :stroke)
    end 50 50
end

for (T, N) in [(X, "X"), (Y, "Y"), (Z, "Z"), (H, "H"), (S, "S"), (Sd, L"S^\dagger"), (T, "T"), (Td, L"T^\dagger")]
    @eval begin
        draw(::$T; kwargs...) = draw($T; kwargs...)
        draw(::Type{$T}; kwargs...) = draw_block($N; kwargs...)
    end
end

function draw(gate::Control)
    c = control(gate)
    t = target(gate)
    r = range(extrema(lanes(gate))...)

    @assert length(t) == 1

    vcat(
        [
            if lane == first(r)
                draw_copy(:top)
            elseif lane ∈ c
                draw_copy(:mid)
            else
                draw_cross()
            end for lane in r if lane < only(t)
        ]...,
        draw(op(gate); top = !any(<(only(t)), c), bottom = !any(>(only(t)), c)),
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
    @drawsvg begin
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

function draw_multiblock_mid()
    @drawsvg begin
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action = :stroke)
        line(Point(25, 0), Point(15, 0), action = :stroke)

        # vertical lines
        line(Point(-15, -25), Point(-15, 25), action = :stroke)
        line(Point(15, -25), Point(15, 25), action = :stroke)
    end 50 50
end

function draw_cross()
    @drawsvg begin
        origin()

        line(Point(-25, 0), Point(25, 0), action = :stroke)
        line(Point(0, -25), Point(0, 25), action = :stroke)
    end 50 50
end

function draw_copy(dir::Symbol)
    @drawsvg begin
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