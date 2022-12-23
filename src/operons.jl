@component function Operon(promoter::Component, genes::Vector{<:Component}; name)
    express!(promoter, genes)
    extend(ReactionSystem(;name=name), promoter.reaction_system)
end
