module GeneticLogicGraph

using ModelingToolkit: @variables, @parameters, @nonamespace
using ModelingToolkit: get_variables, ParentScope, GlobalScope, getvar
using ModelingToolkit: JumpSystem
using Catalyst: @reaction_network, @reaction
using Catalyst: ReactionSystem, Reaction, addreaction!

import ModelingToolkit: equations, states, parameters, get_states, flatten

abstract type Component end
export Component

# methods from ModelingToolkit
equations(x::Component)  = equations(x.reaction_system)
parameters(x::Component) = parameters(x.reaction_system)
states(x::Component)     = states(x.reaction_system)
get_states(x::Component) = get_states(x.reaction_system)
flatten(x::Component)    = flatten(x.reaction_system)

# methods from Base
Base.convert(::Type{ReactionSystem}, x::Component) = x.reaction_system
Base.nameof(x::Component) = nameof(x.reaction_system)
Base.show(io::IO, x::Component) = show(io, x.reaction_system)

function Base.show(io::IO, m::MIME"text/plain", x::Component)
    show(io, m, x.reaction_system)
end

function Base.getproperty(x::Component, sym::Symbol)
    sys = Base.getfield(x, :reaction_system)
    sym == :reaction_system && return sys
    return Base.getproperty(sys, sym)
end

macro component(name, super)
    :(abstract type $(esc(name)) <: $(esc(super)) end)
end
macro component(name, super, network)
    :(begin
          struct $(esc(name)) <: $(esc(super))
              reaction_system::ReactionSystem
              $(esc(name))(x::ReactionSystem) = new(x)
          end
          function $(esc(name))(; name)
              rs = $network
              $(esc(name))(rs)
          end
      end)
end

export equations, states, parameters, flatten

include("products.jl")
export Products
export Monomer
export Dimer

include("promoters.jl")
export Promoters
export SimplePromoter
export TwoStatePromoter
export express!

include("graph.jl")
export Circuit

end # module GeneticLogicGraph
