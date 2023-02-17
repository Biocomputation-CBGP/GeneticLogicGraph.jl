function EmptyCircuit(;name)::ReactionSystem
    @variables t N(t)=0
    @parameters μ
    rx = @reaction μ, ∅ --> N
    ReactionSystem(rx, t, [N], [μ]; name=name, connection_type=Circuit)
end

function Circuit(A::Matrix{Bool}, systems; name)
    rxs = Reaction[]
    syss = ReactionSystem[]
    for i in 1:size(A, 1)
        idxs = findall(identity, A[i, :])
        if length(idxs) > 0
            append!(syss, systems[[i; idxs]])
            for sys in systems[[i; idxs]]
                push!(syss, sys)
                append!(syss, collect(component_args(sys)))
            end
            append!(rxs, connections(systems[i], systems[idxs]...))
        end
    end
    @variables t N(t)=0
    @parameters μ
    push!(rxs, (@reaction μ, ∅ --> N))
    return ReactionSystem(
        rxs, t, [N], [μ];
        name=name, connection_type=(Circuit,), systems=unique(syss)
    )
end

function Circuit(G::DiGraph, systems; name)
    Circuit(Matrix(Graphs.adjacency_matrix(G) .> 0), systems; name=name)
end

function find_N_idx(vars::Vector{T}) where {T<:Term}
    return findfirst(x -> nameof(x.f) == :N, vars)
end

function find_dilutable_idxs(vars::Vector{T}) where {T<:Term}
    return collect(1:length(vars))[isdilutable.(vars)]
end

function make_doubling_callback(model)
    N::Int = find_N_idx(states(model))
    I::Vector{Int} = find_dilutable_idxs(states(model))
    condition = let N=N
        (u, t, integrator) -> u[N] > 0
    end
    affect! = let N=N, I=I
        function (integrator)
            integrator.u[I] .= rand.(Binomial.(@view integrator.u[I]))
            integrator.u[N] = 0
            reset_aggregated_jumps!(integrator)
        end
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end

function make_termination_callback(dt_min)
    condition = let dt=dt_min
        (u, t, integrator) -> (integrator.tstop - integrator.t) < dt_min
    end
    function affect!(integrator)
        integrator.u .= typemax(Int)
        terminate!(integrator)
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end
