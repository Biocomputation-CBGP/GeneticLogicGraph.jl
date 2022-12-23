@component function Promoter(; name, k=1)
    k, = @parameters k=k
    t, promoter = @variables t promoter(t)
    ReactionSystem(Reaction[], t, [promoter], [k]; name=name)
end

express(p::Promoter, x::Vector{<:Component}) = return express(p.reaction_system, x)
function express(p::ReactionSystem, x::Vector{<:Component})
    N = length(x)
    if N > 0
        promoter = @nonamespace p.promoter
        k = @nonamespace p.k
        rnas = GlobalScope.(getproperty.(x, :rna))
        products = [promoter; rnas]
        rxs = [Reaction(k, [promoter], products, [1], fill(1, N+1))]
        p = Promoter(ReactionSystem(rxs, p.iv, states(p), parameters(p); name=nameof(p)))
    end
    return p
end

random_state(x::Promoter) = Dict(x.promoter => rand(0:1))

@component function TwoStatePromoter(a::Promoter, b::Promoter, factor; name, binding=1, unbinding=1)
    @variables t
    @parameters binding=binding unbinding=unbinding
    factor = GlobalScope(factor)
    rxs = [
        (@reaction $(binding), $(a.promoter) + $factor --> $(b.promoter)),
        (@reaction $(unbinding), $(b.promoter) --> $(a.promoter) + $factor),
    ]
    systems = [a.reaction_system, b.reaction_system]
    ReactionSystem(rxs, t, [], [binding, unbinding]; systems=systems, name=name)
end

function express(p::TwoStatePromoter, x::Vector{<:Component})
    a = express(p.systems[1], x).reaction_system
    b = express(p.systems[2], x).reaction_system
    return TwoStatePromoter(
        ReactionSystem(reactions(p), p.iv, get_states(p), p.ps; name=nameof(p), systems=[a,b]))
end

function random_state(x::TwoStatePromoter)
    i = rand(0:1)
    p1 = namespace_expr(x.systems[1].promoter, x.reaction_system, Symbol(:promoters, :₊, nameof(x)))
    p2 = namespace_expr(x.systems[2].promoter, x.reaction_system, Symbol(:promoters, :₊, nameof(x)))
    Dict(p1 => i, p2 => 1 - i)
end



