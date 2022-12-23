@component function Circuit(edges, promoters, products; name)
    N = length(promoters)
    M = falses(N, N)
    for edge ∈ edges
        M[first(edge), last(edge)] = true
    end

    used_promoters = Component[]
    used_products  = falses(N)
    
    for i ∈ 1:N
        outs = M[i, :]
        if any(outs)
            used_products[outs] .= true
            used_products[i] = true
            promoter = express(promoters[i], products[outs])
            push!(used_promoters, promoter)
        end
    end

    @variables t N(t)=0
    @parameters growth=1/45
    rs = @reaction_network $name begin
        $(growth), ∅ --> $N
    end

    products  = convert.(ReactionSystem, products[used_products])
    promoters = convert.(ReactionSystem, used_promoters)
    @named promoters = ReactionSystem(Reaction[], rs.iv; systems=promoters)
    compose(rs, [products; promoters])
end

flatten(x::Circuit) = Circuit(flatten(x.reaction_system))
reactioncomplexes(x::Circuit) = reactioncomplexes(x.reaction_system)
subnetworks(x::Circuit) = subnetworks(x.reaction_system)
promoters(x::Circuit) = x.promoters.systems
function growable(x::Circuit)
    reduce(vcat, [growable(sys) for sys in x.systems], init=typeof(states(x))[])
end
function growable(x::ReactionSystem)
    if nameof(x) != :promoters && !startswith(string(nameof(x)), "constant")
        states(x, states(x))
    else
        typeof(states(x))[]
    end
end

function to_matrix(x::Circuit)
    rs = x.reaction_system
    N = length(states(x))
    involved_in = (substoichmat(rs) .> 0) .| (prodstoichmat(rs) .> 0)
    M = falses(N, N)
    for i ∈ 1:N
        for j ∈ 1:N
            interacts = any((@view involved_in[i, :]) .& (@view involved_in[j, :]))
            M[i, j] = interacts
        end
    end
    M
end

function connected!(from::Int, M::BitMatrix, visited::Vector{Int})
    if from ∉ visited
        N, _ = size(M)
        push!(visited, from)
        nexts = collect(1:N)[@view M[from, :]]
        for next in nexts
            connected!(next, M, visited)
        end
    end
end

function connected(from, x::Circuit)
    indexof(v) = findfirst(isequal(v), states(x))

    idxs = Int[]
    fromidx = indexof(from)
    if !(fromidx === nothing)
        connected!(fromidx, to_matrix(x), idxs)
    end
    states(x)[idxs]
end


function prune(x::Circuit, inputs, outputs)
    x = flatten(x)
    involved = [@nonamespace x.N]

    includesoutput(vs) = any(any(isequal(x).(vs)) for x in outputs)
    for input in inputs
        components = connected(input, x)
        if includesoutput(components)
            for component in components
                if !any(isequal(component), involved) && includesoutput(connected(component, x))
                    push!(involved, component)
                end
            end
        end
    end

    rxs = Reaction[]
    for rx in equations(x)
        species = unique(vcat(rx.products, rx.substrates))
        if length(species ∩ involved) > 0
            push!(rxs, rx)
        end
    end
    rs = Circuit(ReactionSystem(rxs, x.iv; name=nameof(x), defaults=get_defaults(x)))
    if includesoutput(states(rs))
        rs
    else
        nothing
    end
end

function DiscreteProblem(x::Circuit, args...; kwargs...)
    DiscreteProblem(x.reaction_system, args...; kwargs...)
end
function JumpProblem(x::Circuit, args...; kwargs...)
    JumpProblem(x.reaction_system, args...; kwargs...)
end

function doubling_callback(x::Circuit, jumpproblem::JumpProblem)
    indexof(v) = findfirst(isequal(v), states(x))
    idxs = indexof.(growable(x))
    nidx = indexof(@nonamespace x.N)

    condition(u, t, integrator) = u[nidx] > 0
    function affect!(integrator)
        integrator.u[nidx] = 0
        integrator.u[idxs] .= rand.(Binomial.(integrator.u[idxs]))
        reset_aggregated_jumps!(integrator)
    end
    DiscreteCallback(condition, affect!)
end

function sampling(x::Circuit, jumpproblem::JumpProblem, input_func, outputs, N::Int)    
    callback = doubling_callback(x, jumpproblem)
    if N > 1
        problem = EnsembleProblem(jumpproblem)
        solve(
            problem,
            SSAStepper(),
            EnsembleThreads(),
            trajectories=N,
            callback=callback,
            saveat=1,
        )
    else
        solve(
            jumpproblem,
            SSAStepper(),
            callback=callback,
            saveat=1,
        )
    end
end
