module GeneticLogicGraph

using Symbolics
using SymbolicUtils: Symbolic
using ModelingToolkit
using Catalyst
using Distributions
using JumpProcesses
using Graphs


struct Dilute end
Symbolics.option_to_metadata_type(::Val{:dilute}) = Dilute
isdilutable(x)::Bool = getmetadata(x, Dilute, false)
export isdilutable

abstract type Species end

function ConcreteSystemType(::Type{T}) where {T<:Species}
    return ReactionSystem{
        Nothing,
        Catalyst.NetworkProperties{Int, Term{Real, Base.ImmutableDict{DataType, Any}}}
    }
end

abstract type InputSpecies <: Species end
abstract type Monomer <: Species end
abstract type Dimer <: Species end
promote_rule(::Type{InputSpecies}, ::Type{Monomer}) = Species
promote_rule(::Type{InputSpecies}, ::Type{Dimer}) = Species
promote_rule(::Type{Monomer}, ::Type{Dimer}) = Species
component_type(x) = first(x.connection_type)
component_args(x) = x.connection_type[2:end]

export InputSpecies
export Monomer
export Dimer

include("utils.jl")
export prune

include("products.jl")
export randu0

abstract type PromoterRegion end

function ConcreteSystemType(::Type{T}) where {T<:PromoterRegion}
    return ReactionSystem{
        Nothing,
        Catalyst.NetworkProperties{Int, Term{Real, Base.ImmutableDict{DataType, Any}}}
    }
end

abstract type RegulatedPromoter <: PromoterRegion end

export PromoterRegion
export RegulatedPromoter

include("promoters.jl")
export randu0

include("operons.jl")
export connections

abstract type Circuit end
export Circuit

include("graph.jl")
export SingleEdge
export make_doubling_callback
export make_termination_callback

end # module GeneticLogicGraph
