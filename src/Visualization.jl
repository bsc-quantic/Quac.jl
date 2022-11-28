using Luxor
using MathTeXEngine

abstract type VisualizationLayout <: Blueprint end

function compile(::Type{<:VisualizationLayout}, circuit)
    new = Circuit(lanes(circuit))

    for gate in circuit
        # TODO
    end

    error("not implemented yet")
end

function Base.show(io::IO, ::MIME"image/svg+xml", circ::Circuit)
    error("not implemented yet")
end

function draw end
export draw

function draw(gate::AbstractGate; top::Bool=false, bottom::Bool=false)
    n = (length ∘ lanes)(gate)
    if n == 1
        draw_block(; top=bottom, bottom=top)
    else
        a, b = extrema(lanes(gate))
        n = b - a + 1
        vcat(
            draw_multiblock_top(; top=top),
            fill(draw_multiblock_mid(), (n - 2))...,
            draw_multiblock_bottom(; bottom=bottom)
        )
    end
end

draw(::I) = draw(I)
function draw(::Type{I})
    @drawsvg begin
        origin()
        line(Point(-25, 0), Point(25, 0), action=:stroke)
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
            end
            for lane in r if lane < only(t)
        ]...,
        draw(op(gate); top=any(<(only(t)), c), bottom=any(>(only(t)), c)),
        [
            if lane == last(r)
                draw_copy(:bottom)
            elseif lane ∈ c
                draw_copy(:mid)
            else
                draw_cross()
            end
            for lane in r if lane > only(t)
        ]...
    )
end

function draw_block(label=""; top::Bool=false, bottom::Bool=false)
    @drawsvg begin
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action=:stroke)
        line(Point(25, 0), Point(15, 0), action=:stroke)

        # control connectors
        if top
            line(Point(0, 25), Point(0, 15), action=:stroke)
        end
        if bottom
            line(Point(0, -25), Point(0, -15), action=:stroke)
        end

        rect(-15, -15, 30, 30, action=:stroke)

        # label
        fontsize(16)
        text(label, valign=:middle, halign=:center)
    end 50 50
end

function draw_multiblock_mid()
    @drawsvg begin
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action=:stroke)
        line(Point(25, 0), Point(15, 0), action=:stroke)

        # vertical lines
        line(Point(-15, -25), Point(-15, 25), action=:stroke)
        line(Point(15, -25), Point(15, 25), action=:stroke)
    end 50 50
end

function draw_multiblock_bottom()
    @drawsvg begin
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action=:stroke)
        line(Point(25, 0), Point(15, 0), action=:stroke)

        # vertical lines
        line(Point(-15, -25), Point(-15, 15), action=:stroke)
        line(Point(15, -25), Point(15, 15), action=:stroke)

        # terminal
        line(Point(-15, 15), Point(15, 15), action=:stroke)
    end 50 50
end

function draw_multiblock_top()
    @drawsvg begin
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action=:stroke)
        line(Point(25, 0), Point(15, 0), action=:stroke)

        # vertical lines
        line(Point(-15, 25), Point(-15, -15), action=:stroke)
        line(Point(15, 25), Point(15, -15), action=:stroke)

        # terminal
        line(Point(-15, -15), Point(15, -15), action=:stroke)
    end 50 50
end

function draw_cross()
    @drawsvg begin
        origin()

        line(Point(-25, 0), Point(25, 0), action=:stroke)
        line(Point(0, -25), Point(0, 25), action=:stroke)
    end 50 50
end

function draw_copy(dir::Symbol)
    @drawsvg begin
        origin()
        line(Point(-25, 0), Point(25, 0), action=:stroke)

        circle(0, 0, 5, action=:fill)

        (a, b) = if dir == :top
            Point(0, 0), Point(0, 25)
        elseif dir == :bottom
            Point(0, 0), Point(0, -25)
        elseif dir == :mid
            Point(0, 25), Point(0, -25)
        else
            throw(ArgumentError("`dir`=$dir is invalid"))
        end
        line(a, b, action=:stroke)
    end 50 50
end