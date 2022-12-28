import LinearAlgebra

"""
    commutes(A, B)::Bool
"""
function commute end
export commute

commute(A::T, B::T) where {T<:Gate} = true

commute(A::AbstractMatrix, B::AbstractMatrix; kwargs...) = ≈(LinearAlgebra.norm(commutator(A, B)), 0; kwargs...)
function commute(A::AbstractArray{<:Any,N}, B::AbstractArray{<:Any,N}; kwargs...) where {N}
    A = reshape(A, prod(size(A)[1:N÷2]), prod(size(A)[N÷2+1:end]))
    B = reshape(B, prod(size(B)[1:N÷2]), prod(size(B)[N÷2+1:end]))
    commute(A, B; kwargs...)
end

# TODO maybe cache results?
# TODO maybe use `@nospecialize`?
function commute(A::Gate, B::Gate)
    # if acting in disjoint lanes, trivially commute
    isdisjoint(lanes(A), lanes(B)) && return true

    # if acting on same gates, use commutation rules
    lanes(A) == lanes(B) && return commute(typeof(A), typeof(B))

    # else: some overlap on acting lanes
    error("not implemented yet")
end

commute(A::Type{<:Gate}, B::Type{<:Gate}) = commute(arraytype(A)(A), arraytype(B)(B))

commute(::Type{I}, ::Type{<:Gate}) = true
commute(::Type{<:Gate}, ::Type{I}) = true

commute(::Type{X}, ::Type{<:Union{Y,Z}}) = false
commute(::Type{Y}, ::Type{<:Union{X,Z}}) = false
commute(::Type{Z}, ::Type{<:Union{X,Y}}) = false

commute(::Type{H}, ::Type{<:Pauli}) = false
commute(::Type{<:Pauli}, ::Type{H}) = false

commute(::Type{<:Union{X,Rx}}, ::Type{<:Union{X,Rx}}) = true
commute(::Type{<:Union{Y,Ry}}, ::Type{<:Union{Y,Ry}}) = true
commute(::Type{<:Phase}, ::Type{<:Phase}) = true

"""
    commutator(A,B)

Computes the commutator of `A` with `B`.
"""
function commutator(A, B)
    A * B - B * A
end

"""
    anticommutator(A,B)

Computes the anticommutator of `A` with `B`.
"""
function anticommutator(A, B)
    A * B + B * A
end