function PromoterRegion(transcription; name)
    @variables t
    @variables promoter(t)      [description="RNAP recruiting region of DNA", dilute=false]
    @parameters λ=transcription [description="Transcription rate from this promoter region"]
    opts = Dict(:name => name, :connection_type => (PromoterRegion, ))
    return ReactionSystem(Reaction[], t, [promoter], [λ]; opts...)
end

function RegulatedPromoter(bound::ReactionSystem, unbound::ReactionSystem, ligand::ReactionSystem, binding, unbinding; name)
    Rs = GlobalScope.(states(ligand, ModelingToolkit.outputs(ligand)))
    @variables t
    @parameters k₁=binding   [description="Binding rate of the ligand to the promoter region"]
    @parameters k₀=unbinding [description="Unbinding rate of the ligand from the promoter region"]
    bindings   = [Reaction(k₁, [unbound.promoter, R], [bound.promoter]) for R in Rs]
    unbindings = [Reaction(k₀, [bound.promoter], [R, unbound.promoter]) for R in Rs]
    opts = Dict(
        :name => name,
        :systems => [unbound, bound],
        :connection_type => (RegulatedPromoter, ligand)
    )
    return ReactionSystem([bindings; unbindings], t, [], [k₀, k₁]; opts...)
end

function RegulatedPromoter(bound::Real, unbound::Real, ligand::ReactionSystem, binding, unbinding; name)
    @named bound   = PromoterRegion(bound)
    @named unbound = PromoterRegion(unbound)
    RegulatedPromoter(bound, unbound, ligand, binding, unbinding; name=name)
end

randu0(::Type{PromoterRegion}, x::ReactionSystem) = Dict(x.promoter => 1)
function randu0(::Type{RegulatedPromoter}, x::ReactionSystem)
    i = rand(0:1)
    Dict(x.bound.promoter => i, x.unbound.promoter => 1 - i)
end
