function InputSpecies(production, degradation; name)
    @parameters λ=production  [description="Rate of production for this species"]
    @parameters α=degradation [description="Rate of degradation for this species"]

    @variables t
    @variables species(t) [
        description="The abundance of the species",
        dilute=true,
        output=true
    ]
    rxs = [
        Reaction(λ, nothing, [species]),
        Reaction(α, [species], nothing),
    ]
    opts = Dict(:name => name, :connection_type => (InputSpecies, ))
    return ReactionSystem(rxs, t, [species], [λ, α]; opts...)
end

function Monomer(translation, rna_degradation, degradation; name)
    @parameters λ=translation      [description="Rate of protein synthesis (translation) from mRNA"]
    @parameters α₁=rna_degradation [description="Rate of RNA degradation"]
    @parameters α₂=degradation     [description="Rate of protein degradation"]
    
    @variables t
    @variables rna(t) [
        description="Abundance of mRNA for the monomer",
        dilute=true,
        input=true
    ]
    @variables monomer(t) [
        description="Abundance of monomers",
        dilute=true,
        output=true
    ]
    rxs = [
        Reaction(λ, [rna], [rna, monomer]),
        Reaction(α₁, [rna], nothing),
        Reaction(α₂, [monomer], nothing),
    ]
    opts = Dict(:name => name, :connection_type => (Monomer, ))
    return ReactionSystem(rxs, t, [rna, monomer], [λ, α₁, α₂]; opts...)
end

function Dimer(translation, rna_degradation, degradation, binding, unbinding; name)
    @parameters λ=translation      [description="Rate of protein synthesis (translation) from mRNA"]
    @parameters α₁=rna_degradation [description="Rate of RNA degradation"]
    @parameters α₂=degradation     [description="Rate of protein degradation"]
    @parameters k₋₁=unbinding      [description="Rate of dimer unbinding"]
    @parameters k₁=binding         [description="Rate of monomer dimerisation"]
    
    @variables t
    @variables rna(t) [
        description="Abundance of mRNA for the monomer",
        dilute=true,
        input=true
    ]
    @variables monomer(t) [
        description="Abundance of monomers",
        dilute=true,
    ]
    @variables dimer(t) [
        description="Abundance of dimers",
        dilute=true,
        output=true
    ]
    rxs = [
        Reaction(λ, [rna], [rna, monomer]),
        Reaction(α₁, [rna], []),
        Reaction(α₂, [monomer], []),
        Reaction(α₂, [dimer], []),
        Reaction(k₁, [monomer], [dimer], [2], [1]),
        Reaction(k₋₁, [dimer], [monomer], [1], [2]),
    ]
    T =  ConcreteSystemType(Dimer)
    opts = Dict(:name => name, :connection_type => (Dimer, ))
    return ReactionSystem(rxs, t, [rna, monomer, dimer], [λ, α₁, α₂, k₋₁, k₁]; opts...)
end

randu0(x::ReactionSystem) = randu0(component_type(x), x)
randu0(::Type{InputSpecies}, x) = Dict(x.species => rand(0:16))
randu0(::Type{Monomer}, x) = Dict(x.rna => rand(0:16), x.monomer => rand(0:16))
randu0(::Type{Dimer}, x) = merge(randu0(Monomer, x), Dict(x.dimer => rand(0:4)))

output(x::ReactionSystem) = output(component_type(x), x)
output(::Type{InputSpecies}, x::ReactionSystem) = x.species
output(::Type{Monomer}, x::ReactionSystem) = x.monomer
output(::Type{Dimer}, x::ReactionSystem) = x.dimer
