struct Circuit <: Component
    reaction_system::ReactionSystem
    Circuit(x::ReactionSystem) = new(x)
end
function Circuit(edges, mapping; name)
    N = length(mapping)
    M = falses(N, N)
    for edge ∈ edges
        M[first(edge), last(edge)] = true
    end
    products = convert(ReactionSystem, Products(unique(last.(mapping))))
    all_promoters = first.(mapping)

    for i ∈ 1:N
        tfs = last.(mapping[M[i, :]])
        rnas = []
        for tf in tfs
            sys = getvar(products, nameof(tf); namespace=true)
            push!(rnas, getvar(sys, :rna; namespace=true))
        end
        all_promoters[i] = express!(all_promoters[i], rnas...)
    end


    promoters = convert(ReactionSystem, Promoters(all_promoters))
    rs = ReactionSystem([@reaction growth, ∅ --> N]; name=name, systems=[products, promoters])
    Circuit(rs)
end

