module QuacMakieExt

if isdefined(Base, :get_extension)
    using Quac
else
    using ..Quac
end

using Makie

@recipe(Blueprint, circuit) do scene
    labels_theme = default_theme(scene, Makie.text)
    Attributes(; default_theme(scene)..., lanecolor = :black, lanewidth = 2)
end

Makie.plottype(::Circuit) = Blueprint

function Makie.plot!(bp::Blueprint)
    circuit = bp[:circuit]

    n = lanes(circuit[])
    m = max(1, length(moments(circuit[])))

    hlines!(bp, 1:n, xmax = m, color = bp[:lanecolor], linewidth = bp[:lanewidth])

    for gate in circuit[]
        poly!(bp, Rect())
    end
end

function draw! end

end
