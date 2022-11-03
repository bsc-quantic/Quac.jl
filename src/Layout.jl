abstract type Layout end

function transform! end

function transform(::Type{L}, circ::Circuit) where {L<:Layout}
    new_circ = deepcopy(circ)
    transform!(L, circ)

    return circ
end

transform!(::Type{<:Layout}, _::Circuit) = error("not yet implemented")


struct NonOverlappingLayout <: Layout end
