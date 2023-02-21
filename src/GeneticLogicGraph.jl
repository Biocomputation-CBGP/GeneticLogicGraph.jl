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

struct mRNA end
Symbolics.option_to_metadata_type(::Val{:mrna}) = mRNA
ismrna(x)::Bool = getmetadata(x, mRNA, false)
export ismrna

struct Protein end
Symbolics.option_to_metadata_type(::Val{:protein}) = Protein
isprotein(x)::Bool = getmetadata(x, Protein, false)
export isprotein

abstract type Component end
function ConcreteSystemType(::Type{T}) where {T<:Component}
    return ReactionSystem{
        Nothing,
        Catalyst.NetworkProperties{Int, Term{Real, Base.ImmutableDict{DataType, Any}}}
    }
end

abstract type Species <: Component end
abstract type InputSpecies <: Species end
abstract type ConstantSpecies <: Species end
abstract type Monomer <: Species end
abstract type Dimer <: Species end
promote_rule(::Type{InputSpecies}, ::Type{Monomer}) = Species
promote_rule(::Type{InputSpecies}, ::Type{Dimer}) = Species
promote_rule(::Type{ConstantSpecies}, ::Type{Monomer}) = Species
promote_rule(::Type{ConstantSpecies}, ::Type{Dimer}) = Species
promote_rule(::Type{Monomer}, ::Type{Dimer}) = Species
component_type(x) = first(x.connection_type)
component_args(x) = x.connection_type[2:end]

translation_reactions(x::ReactionSystem, args...) = translation_reactions(component_type(x), x, args...)
translation_reactions(::Type{<:Component}, args...) = Reaction[]

function transcription_reactions(x::T, args::Vararg{T}) where {T<:ReactionSystem}
    T1 = component_type(x)
    T2 = promote_type(component_type.(args)...)
    return transcription_reactions(T1, T2, x, args...)
end
transcription_reactions(::Type{<:Component}, args...) = Reaction[]

mrna_degradation_reactions(x::ReactionSystem, args...) = mrna_degradation_reactions(component_type(x), x, args...)
mrna_degradation_reactions(::Type{<:Component}, args...) = Reaction[]

export InputSpecies
export ConstantSpecies
export Monomer
export Dimer

include("utils.jl")
export prune

include("products.jl")
export randu0
export translation_reactions
export mrna_degradation_reactions

abstract type PromoterRegion <: Component end
abstract type RegulatedPromoter <: PromoterRegion end

export PromoterRegion
export RegulatedPromoter

include("promoters.jl")
export randu0
export transcription_reactions

abstract type Circuit end
export Circuit

include("graph.jl")
export make_doubling_callback
export make_termination_callback

end # module GeneticLogicGraph
