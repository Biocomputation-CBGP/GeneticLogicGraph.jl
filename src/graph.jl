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
    rs = ReactionSystem(rxs, t, [N], [μ]; name=name, connection_type=(Circuit,))
    return compose(rs, unique(syss))
end

function Circuit(G::DiGraph, systems; name)
    Circuit(Matrix(Graphs.adjacency_matrix(G) .> 0), systems; name=name)
end

function make_doubling_callback(model)
    idxs = findall(isdilutable, states(model))
    Nidx = findfirst(isequal(@nonamespace model.N), states(model))

    function affect!(integrator)
        integrator.u[idxs] .= rand.(Binomial.(@view integrator.u[idxs]))
        integrator.u[Nidx] = 0
        reset_aggregated_jumps!(integrator)
    end
    condition(u, t, integrator) = u[Nidx] == 1
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end

function make_termination_callback(dt_min)
    condition(u, t, integrator) = (integrator.tstop - integrator.t) < dt_min
    function affect!(integrator)
        @debug "The simulation was terminated because dt < $(dt_min)"
        integrator.u .= -1
        terminate!(integrator)
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end
