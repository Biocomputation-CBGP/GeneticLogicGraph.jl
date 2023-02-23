"""
    InputSpecies(production, degradation; name)

Construct a model of a chemical species which is created from nothing.

InputSpecies can be used as "sources" in larger models. They have no
transcription or translation reactions, only production and
degradation.

# Example
```julia-repl
julia> using ModelingToolkit

julia> @named IPTG = InputSpecies(1, 0.1)
```
"""
function InputSpecies(production; name)
    @parameters λ=production  [description="Rate of production for this species"]
    @variables t
    @variables species(t) [
        description="The abundance of the species",
        dilute=true,
        output=true
    ]
    rxs = [
        Reaction(λ, nothing, [species]),
    ]
    opts = Dict(:name => name, :connection_type => (InputSpecies, ))
    return ReactionSystem(rxs, t, [species], [λ]; opts...)
end

"""
    ConstantSpecies(level; name)

Construct a model of a chemical species which is has constant abundance.

ConstantSpecies can be used as "sources" in larger models. They are
not produced or degraded or diluted.

# Example
```julia-repl
julia> using ModelingToolkit

julia> @named IPTG = ConstantSpecies(1)
```
"""
function ConstantSpecies(level; name)
    @variables t
    @variables species(t)=level [
        description="The abundance of the species",
        dilute=false,
        output=true
    ]
    opts = Dict(:name => name, :connection_type => (ConstantSpecies, ))
    return ReactionSystem(Reaction[], t, [species], []; opts...)
end

function Monomer(translation; name)
    @parameters λ=translation [description="Rate of protein synthesis (translation) from mRNA"]
    @variables t
    @variables rna(t) [
        description="Abundance of mRNA for the monomer",
        dilute=true,
        mrna=true,
        input=true
    ]
    @variables monomer(t) [
        description="Abundance of monomers",
        dilute=true,
        protein=true,
        output=true
    ]    
    opts = Dict(:name => name, :connection_type => (Monomer, ))
    return ReactionSystem(Reaction[], t, [rna, monomer], [λ]; opts...)
end

function Dimer(translation, binding, unbinding; name)
    @parameters λ=translation      [description="Rate of protein synthesis (translation) from mRNA"]
    @parameters k₋₁=unbinding      [description="Rate of dimer unbinding"]
    @parameters k₁=binding         [description="Rate of monomer dimerisation"]
    
    @variables t
    @variables rna(t) [
        description="Abundance of mRNA for the monomer",
        mrna=true,
        dilute=true,
        input=true
    ]
    @variables monomer(t) [
        description="Abundance of monomers",
        protein=true,
        dilute=true,
    ]
    @variables dimer(t) [
        description="Abundance of dimers",
        protein=true,
        dilute=true,
        output=true
    ]
    rxs = [
        Reaction(k₁, [monomer], [dimer], [2], [1]),
        Reaction(k₋₁, [dimer], [monomer], [1], [2]),
    ]
    opts = Dict(:name => name, :connection_type => (Dimer, ))
    return ReactionSystem(rxs, t, [rna, monomer, dimer], [λ, k₋₁, k₁]; opts...)
end

function _translation_reactions(x::ReactionSystem)
    return [Reaction(x.λ, [x.rna], [x.rna, x.monomer])]
end

function _translation_reactions(x::ReactionSystem, ribosome, r₁, r₋₁)
    s = Symbol((@nonamespace x.rna).val.f.name, "_ribosome_complex")
    vs = @variables t $s(t)=0 [
        description="The ribsosome-rna complex",
        dilute=true
    ]
    complex = vs[2]
    addspecies!(x, complex)
    return [
        Reaction(r₁, [x.rna, ribosome], [getproperty(x, s)]),
        Reaction(r₋₁, [getproperty(x, s)], [x.rna, ribosome]),
        Reaction(x.λ, [getproperty(x, s)], [x.rna, ribosome, x.monomer]),
    ]
end

function _mrna_degradation_reactions(x::ReactionSystem, α)
    return [Reaction(α, [x.rna], nothing)]
end

translation_reactions(::Type{Monomer}, x::ReactionSystem, args...) = _translation_reactions(x, args...)
translation_reactions(::Type{Dimer}, x::ReactionSystem, args...) = _translation_reactions(x, args...)
mrna_degradation_reactions(::Type{Monomer}, x::ReactionSystem, args...) = _mrna_degradation_reactions(x, args...)
mrna_degradation_reactions(::Type{Dimer}, x::ReactionSystem, args...) = _mrna_degradation_reactions(x, args...)

randu0(x::ReactionSystem) = randu0(component_type(x), x)
randu0(::Type{ConstantSpecies}, x) = Dict(x.species => rand(0:16))
randu0(::Type{InputSpecies}, x) = Dict(x.species => rand(0:16))
randu0(::Type{Monomer}, x) = Dict(x.rna => rand(0:16), x.monomer => rand(0:16))
randu0(::Type{Dimer}, x) = merge(randu0(Monomer, x), Dict(x.dimer => rand(0:4)))

output(x::ReactionSystem) = output(component_type(x), x)
output(::Type{InputSpecies}, x::ReactionSystem) = x.species
output(::Type{Monomer}, x::ReactionSystem) = x.monomer
output(::Type{Dimer}, x::ReactionSystem) = x.dimer
