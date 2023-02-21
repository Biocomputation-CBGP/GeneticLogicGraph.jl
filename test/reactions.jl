using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using Test

const GLG = GeneticLogicGraph


@testset "Connections between promoters" begin
    @named input = InputSpecies(1.0, 1.0)
    @named monomer = Monomer(1.0)
    @named dimer = Dimer(1.0, 1.0, 1.0)

    @named simplepromoter = PromoterRegion(1.0)
    @named promoter1 = RegulatedPromoter(1.0, 1.0, input, 1.0, 1.0)
    @named promoter2 = RegulatedPromoter(1.0, 1.0, monomer, 1.0, 1.0)
    @named promoter3 = RegulatedPromoter(1.0, 1.0, dimer, 1.0, 1.0)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, monomer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, monomer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(promoter1, promoter2)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, dimer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, dimer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(promoter1, promoter3)
    
    rxs = [
        Reaction(simplepromoter.λ, [simplepromoter.promoter], [simplepromoter.promoter, monomer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(simplepromoter, promoter2)

    rxs = [
        Reaction(simplepromoter.λ, [simplepromoter.promoter], [simplepromoter.promoter, dimer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(simplepromoter, promoter3)

    rxs = [
        Reaction(promoter2.bound.λ, [promoter2.bound.promoter], [promoter2.bound.promoter, monomer.rna]),
        Reaction(promoter2.unbound.λ, [promoter2.unbound.promoter], [promoter2.unbound.promoter, monomer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(promoter2, promoter2)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, monomer.rna, dimer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, monomer.rna, dimer.rna]),
    ]
    @test rxs == GLG.transcription_reactions(promoter1, promoter2, promoter3)
end

@testset "Connections between other components (empty reaction lists)" begin
    @named input = InputSpecies(1.0, 1.0)
    @named monomer = Monomer(1.0)
    @named dimer = Dimer(1.0, 1.0, 1.0)

    @test Reaction[] == GLG.transcription_reactions(input, monomer)
    @test Reaction[] == GLG.transcription_reactions(monomer, input)
    @test Reaction[] == GLG.transcription_reactions(dimer, monomer)
    @test Reaction[] == GLG.transcription_reactions(monomer, dimer)
    @test Reaction[] == GLG.transcription_reactions(input, input)
    @test Reaction[] == GLG.transcription_reactions(monomer, monomer)
    @test Reaction[] == GLG.transcription_reactions(dimer, dimer)
end

@testset "Translation for products" begin
    @named monomer = Monomer(1.0)
    @named dimer = Dimer(1.0, 1.0, 1.0)

    expected = [Reaction(monomer.λ, [monomer.rna], [monomer.rna, monomer.monomer])]
    @test expected == GLG.translation_reactions(monomer)

    expected = [Reaction(dimer.λ, [dimer.rna], [dimer.rna, dimer.monomer])]
    @test expected == GLG.translation_reactions(dimer)

    @variables t ribosome(t) monomer₊rna_ribosome_complex(t) dimer₊rna_ribosome_complex(t)
    @parameters r₁ r₋₁
    expected = [
        Reaction(r₁, [monomer.rna, ribosome], [monomer₊rna_ribosome_complex]),
        Reaction(r₋₁, [monomer₊rna_ribosome_complex], [monomer.rna, ribosome]),
        Reaction(monomer.λ, [monomer₊rna_ribosome_complex], [monomer.rna, ribosome, monomer.monomer])
    ]
    @test expected == GLG.translation_reactions(monomer, ribosome, r₁, r₋₁)
end

@testset "Translation in other components (empty reaction lists)" begin
    @named input = InputSpecies(1.0, 1.0)
    @named constant = ConstantSpecies(1)
    @named promoter = PromoterRegion(1)
    @named regulatedpromoter = RegulatedPromoter(1, 1, input, 1, 1)

    @test Reaction[] == GLG.translation_reactions(input)
    @test Reaction[] == GLG.translation_reactions(constant)
    @test Reaction[] == GLG.translation_reactions(promoter)
    @test Reaction[] == GLG.translation_reactions(regulatedpromoter)
end

@testset "Degradation of mRNA for products" begin
    @named monomer = Monomer(1.0)
    @named dimer = Dimer(1.0, 1.0, 1.0)
    @parameters α

    expected = [Reaction(α, [monomer.rna], nothing)]
    @test expected == GLG.mrna_degradation_reactions(monomer, α)

    expected = [Reaction(α, [dimer.rna], nothing)]
    @test expected == GLG.mrna_degradation_reactions(dimer, α)
end
