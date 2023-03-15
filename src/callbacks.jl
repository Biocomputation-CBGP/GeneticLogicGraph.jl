struct MaxAbundanceError <: Exception
    idx::Int
    abundance::Int
end

function Base.showerror(io::IO, e::MaxAbundanceError)
    msg = "Simulation reached max abundance $(e.abundance) for species $(e.idx)"
    print(io, msg)
    return nothing
end

struct MinTimestepError{T<:Real} <: Exception
    dt::T
    dtmin::T
end

function Base.showerror(io::IO, e::MinTimestepError)
    msg = (
        "Simulation wants to use timestep $(e.dt), "
        * "which is smaller than the minimum timestep $(e.dtmin)"
    )
    print(io, msg)
    return nothing
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

function make_mintimestep_callback(mindt)
    condition = let dt=mindt
        (u, t, integrator) -> integrator.tstop - t < dt
    end
    affect! = let dt=mindt
        (integrator) -> throw(MinTimestepError(integrator.tstop - integrator.t, dt))
    end
    return DiscreteCallback(condition, affect!, save_positions=(false, false))
end
