abstract type Promoter <: Component end

struct Promoters <: Component
    reaction_system::ReactionSystem
    Promoters(x::ReactionSystem) = new(x)
end
function Promoters(promoters::Vector{<:Promoter})
    @variables t
    promoters = [promoter.reaction_system for promoter in promoters]
    Promoters(ReactionSystem(Reaction[], t; systems=promoters, name=:promoters))
end
function Promoters(promoters::Vararg{<:Promoter})
    @variables t
    promoters = [promoter.reaction_system for promoter in promoters]
    Promoters(ReactionSystem(Reaction[], t; systems=promoters, name=:promoters))
end

struct SimplePromoter <: Promoter
    reaction_system::ReactionSystem
    SimplePromoter(x::ReactionSystem) = new(x)
end
function SimplePromoter(; name)
    k, = @parameters k
    t, promoter = @variables t promoter(t)
    SimplePromoter(ReactionSystem(Reaction[], t, [promoter], [k]; name=name))
end

function express!(promoter::SimplePromoter, rnas...)
    express!(promoter.reaction_system, rnas...)
    promoter
end

function express!(promoter::ReactionSystem, rnas...)
    if length(rnas) > 0
        rnas = GlobalScope.(rnas)
        transcription = Reaction(
            (@nonamespace promoter.k),
            [(@nonamespace promoter.promoter)],
            [(@nonamespace promoter.promoter), rnas...],
            [1], ones(Int, 1 + length(rnas)),
        )
        addreaction!(promoter, transcription)
    end
    promoter
end

struct TwoStatePromoter <: Promoter
    reaction_system::ReactionSystem
    TwoStatePromoter(x::ReactionSystem) = new(x)
end
function TwoStatePromoter(p::Promoter, q::Promoter, tf; name)
    ps = @parameters k₁ k₋₁
    @variables t
    tf = GlobalScope(tf)
    rxs = [
        Reaction(k₁, [p.promoter, tf], [q.promoter]),
        Reaction(k₋₁, [q.promoter], [p.promoter, tf]),
    ]
    subsystems = [p.reaction_system, q.reaction_system]
    TwoStatePromoter(ReactionSystem(rxs, t, [], ps; systems=subsystems, name=name))
end

function express!(promoter::TwoStatePromoter, rnas...)
    express!(promoter.systems[1], rnas...)
    express!(promoter.systems[2], rnas...)
    promoter
end



