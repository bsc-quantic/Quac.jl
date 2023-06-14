using Luxor
using MathTeXEngine
using LaTeXStrings

texname(::Type{Op}) where {Op<:Operator} = LaTeXString(String(nameof(Op)))

texname(::Type{Sd}) = L"S^\dagger"
texname(::Type{Td}) = L"T^\dagger"

texname(::Type{Rx}) = L"R_X"
texname(::Type{Ry}) = L"R_Y"
texname(::Type{Rz}) = L"R_Z"

texname(::Type{Hz}) = L"H_Z"
texname(::Type{FSim}) = L"F_S"
function draw end
export draw

function draw(circuit::Circuit; kwargs...)
    n = lanes(circuit)

    if isempty(circuit)
        return vcat([draw(Gate{I}(lane); kwargs...) for lane in 1:n]...)
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

        mapreduce(x -> draw(x; kwargs...), vcat, moment)
    end
end

function Base.show(io::IO, ::MIME"image/svg+xml", circuit::Circuit)
    _ = draw(circuit)
    print(io, svgstring())
end

function draw(gate::Gate{Op,1,P}; kwargs...) where {Op,P}
    draw_block(; top = false, bottom = false)
end

function draw(gate::Gate{Op,N,P}; kwargs...) where {Op,N,P}
    a, b = extrema(lanes(gate))
    n = b - a + 1
    vcat(
        draw_block(; top = true, bottom = false, kwargs...),
        fill(draw_multiblock_mid(; kwargs...), (n - 2))...,
        draw_block(; top = false, bottom = true, kwargs...),
    )
end

function draw(::Gate{I,1,NamedTuple{(),Tuple{}}}; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        origin()
        line(Point(-25, 0), Point(25, 0), action = :stroke)
    end 50 50
end

for Op in [X, Y, Z, H, S, Sd, T, Td, Rx, Ry, Rz, Hz, FSim]
    @eval draw(::Gate{$Op,1,P}; kwargs...) where {P} = draw_block(texname($Op); kwargs...)
end

function draw(gate::Gate{<:Control}; kwargs...)
    c = control(gate)
    t = target(gate)
    r = range(extrema(lanes(gate))...)

    # TODO Control{Swap}
    @assert length(t) == 1

    vcat(
        [
            if lane == first(r)
                draw_copy(:top; kwargs...)
            elseif lane ∈ c
                draw_copy(:mid; kwargs...)
            else
                draw_cross(; kwargs...)
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
                draw_copy(:bottom; kwargs...)
            elseif lane ∈ c
                draw_copy(:mid; kwargs...)
            else
                draw_cross(; kwargs...)
            end for lane in r if lane > only(t)
        ]...,
    )
end

function draw_block(label = ""; top::Bool = false, bottom::Bool = false, background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
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

function draw_multiblock_mid(; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        origin()

        # lane wire
        line(Point(-25, 0), Point(-15, 0), action = :stroke)
        line(Point(25, 0), Point(15, 0), action = :stroke)

        # vertical lines
        line(Point(-25, 0), Point(25, 0), action = :stroke)
        line(Point(0, -25), Point(0, 25), action = :stroke)
    end 50 50
end

function draw_cross(; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
        origin()

        line(Point(-25, 0), Point(25, 0), action = :stroke)
        line(Point(0, -25), Point(0, 25), action = :stroke)
    end 50 50
end

function draw_copy(dir::Symbol; background = nothing)
    @drawsvg begin
        (background !== nothing) && Luxor.background(background)
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