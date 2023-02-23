function transcription_reactions(A::Matrix{Bool}, systems)
    rxs = Reaction[]
    for i in 1:size(A, 1)
        idxs = findall(identity, A[i, :])
        if length(idxs) > 0
            append!(rxs, transcription_reactions(systems[i], systems[idxs]...))
        end
    end
    return rxs
end

function used_systems(A::Matrix{Bool}, systems)
    syss = ReactionSystem[]
    for i in 1:size(A, 1)
        idxs = findall(identity, A[i, :])
        if length(idxs) > 0
            for sys in systems[[i; idxs]]
                push!(syss, sys)
                append!(syss, collect(component_args(sys)))
            end
        end
    end
    return unique(syss)
end

function _FiniteRibosomeCircuit(A::Matrix{Bool}, systems; name)
    syss = used_systems(A, systems)
    @named Ribosome = Monomer(0.1)
    @named pRibosome = PromoterRegion(0.5)
    append!(syss, [Ribosome, pRibosome])

    @variables t N(t)=0
    @parameters μ α r₁ r₋₁

    rxs = [
        transcription_reactions(A, systems);
        transcription_reactions(pRibosome, Ribosome);
        reduce(vcat, mrna_degradation_reactions(sys, α) for sys in syss);
        reduce(vcat, translation_reactions(sys, Ribosome.monomer, r₁, r₋₁) for sys in syss);
        [Reaction(μ, [Ribosome.monomer], [Ribosome.monomer, N])];
    ]

    return ReactionSystem(
        rxs, t, [N], [μ, α, r₁, r₋₁];
        name=name, connection_type=(Circuit,), systems=syss
    )
end

function _Circuit(A::Matrix{Bool}, systems; name)
    syss = used_systems(A, systems)
    @variables t N(t)=0
    @parameters μ α

    rxs = [
        transcription_reactions(A, systems);
        reduce(vcat, mrna_degradation_reactions(sys, α) for sys in syss);
        reduce(vcat, translation_reactions(sys) for sys in syss);
        [Reaction(μ, nothing, [N])]
    ]

    return ReactionSystem(rxs, t, [N], [μ, α]; name=name, connection_type=(Circuit,), systems=syss)
end

function Circuit(A::Matrix{Bool}, systems; name, finite_resource=false)
    if finite_resource
        _FiniteRibosomeCircuit(A, systems; name=name)
    else
        _Circuit(A, systems; name=name)
    end
end

function Circuit(G::DiGraph, systems; name, finite_resource=false)
    Circuit(Matrix(Graphs.adjacency_matrix(G) .> 0), systems; name=name, finite_resource=finite_resource)
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

function make_termination_callback(time_max)
    start_time::Float64 = 0.0
    condition = let time_max=time_max
        function (u, t, integrator)
            if start_time == 0.0
                start_time = time()
            end
            if time() - start_time >= time_max
                start_time = 0.0
                return true
            end
            return false
        end
    end
    function affect!(integrator)
        integrator.u .= 987654321
        @warn "Terminating early"
        terminate!(integrator)
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end
