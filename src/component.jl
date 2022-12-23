abstract type Component end

equations(x::Component)                  = equations(x.reaction_system)
parameters(x::Component)                 = parameters(x.reaction_system)
states(x::Component, args...)            = states(x.reaction_system, args...)
get_states(x::Component)                 = get_states(x.reaction_system)
get_defaults(x::Component)               = get_defaults(x.reaction_system)
reactions(x::Component)                  = reactions(x.reaction_system)
addreaction!(x::Component, rx::Reaction) = addreaction!(x.reaction_system, rx)

random_state(x::Component) = Dict()

Base.convert(::Type{ReactionSystem}, x::Component) = x.reaction_system
function Base.convert(::Type{ReactionSystem}, x::Vector{<:Component})
    convert.(ReactionSystem, x)
end
Base.nameof(x::Component) = nameof(x.reaction_system)
Base.show(io::IO, x::Component) = show(io, x.reaction_system)
function Base.show(io::IO, M::MIME"text/plain", x::Component)
    show(io, M, x.reaction_system)
end

function Base.getproperty(x::Component, sym::Symbol; namespace=true)
    sys = Base.getfield(x, :reaction_system)
    sym == :reaction_system && return sys
    return Base.getproperty(sys, sym; namespace=namespace)
end



macro component(expr)
    head = expr.args[1]
    if head.head == :call
        name = head.args[1]
    elseif head.head == :where
        name = head.args[1].args[1]
    end
    Expr(
        :block, 
        :(struct $(name) <: Component
              reaction_system::ReactionSystem
              $(name)(x::ReactionSystem) = new(x)
          end),
        Expr(:function, esc(expr.args[1]), Expr(:call, esc(name), esc(expr.args[2])))
    )
end
