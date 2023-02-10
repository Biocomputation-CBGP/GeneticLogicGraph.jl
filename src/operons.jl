function connections(x::T, args::Vararg{T}) where {T<:ReactionSystem}
    first_type = component_type(x)
    second_type = promote_type(component_type.(args)...)
    connections(first_type, second_type, x, args...)
end

function connections(::Type{PromoterRegion}, ::Type{<:Species}, x, args...)
    rnas = reduce(vcat, states(y, ModelingToolkit.inputs(y)) for y in args)
    [Reaction(x.λ, [x.promoter], [x.promoter; rnas], [1], ones(Int, length(rnas) + 1))]
end

function connections(::Type{PromoterRegion}, ::Type{RegulatedPromoter}, x, args...)
    connections(x, reduce(vcat, collect(component_args(y)) for y in args)...)
end

function connections(::Type{RegulatedPromoter}, ::Type{RegulatedPromoter}, x, args...)
    vcat(connections(x.bound, args...), connections(x.unbound, args...))
end

function connections(::Type{RegulatedPromoter}, ::Type{<:Species}, x, args...)
    vcat(connections(x.bound, args...), connections(x.unbound, args...))
end
