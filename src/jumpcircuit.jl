function transcription(A::Matrix{Bool}, systems)
    rxs = Reaction[]
    for i in 1:size(A, 1)
        idxs = findall(identity, A[i, :])
        if length(idxs) > 0
            append!(rxs, transcription(systems[i], systems[idxs]...))
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

function Circuit(G::DiGraph, systems; name)
    A = Matrix(Graphs.adjacency_matrix(G) .> 0)
    syss = used_systems(A, systems)
    @variables t N(t)=0
    @parameters μ α

    rxs = [
        transcription(A, systems);
        reduce(vcat, mrna_degradation(sys, α) for sys in syss);
        reduce(vcat, translation(sys) for sys in syss);
        [Reaction(μ, nothing, [N])]
    ]

    return ReactionSystem(
        rxs, t, [N], [μ, α];
        name=name, connection_type=(JumpCircuit,), systems=syss
    )
end

function make_doubling_callback(model)
    @variables t N(t)
    Ni::Int = findfirst(isequal(N), states(model))
    I::Vector{Int} = findall(isdilutable, states(model))
    condition = let N=Ni
        (u, t, integrator) -> u[N] > 0
    end
    affect! = let N=Ni, I=I
        function (integrator)
            integrator.u[I] .= rand.(Binomial.(@view integrator.u[I]))
            integrator.u[N] = 0
            reset_aggregated_jumps!(integrator)
        end
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end

struct MaxAbundanceError <: Exception
    idx::Int
    abundance::Int
end

function Base.showerror(io::IO, e::MaxAbundanceError)
    print(io, "Simulation reached max abundance $(e.abundance) for $(e.idx)")
end

function make_species_limit_callback(max_abundance)
    condition = let max_abundance=max_abundance
        function (u, t, integrator)
            return any(u .> max_abundance)
        end
    end
    function affect!(integrator)
        x, i = findmax(integrator.u)
        throw(MaxAbundanceError(i, x))
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end

function randu0(::Type{JumpCircuit}, x::ReactionSystem)
    return reduce(merge, randu0.(x.systems), init=Dict())
end

function zerou0(::Type{JumpCircuit}, x::ReactionSystem)
    return reduce(merge, zerou0.(x.systems), init=Dict())
end
