using GeneticLogicGraph
using ModelingToolkit
using Catalyst
using Test

const GLG = GeneticLogicGraph

@testset "Promoter region definition" begin
    @named model = PromoterRegion(1.0)

    @test length(equations(model)) == 0
    @test length(states(model)) == 1
    @test length(parameters(model)) == 1

    @variables t promoter(t)
    @test isequal(t, ModelingToolkit.get_iv(model))
    @test any(isequal(promoter), states(model))

    @parameters λ
    @test any(isequal(λ), parameters(model))
    @test ModelingToolkit.get_defaults(model)[λ] == 1.0

    @test GLG.component_type(model) == GLG.PromoterRegion
end

@testset "Regulated Promoter definition" begin
    @named ligand  = InputSpecies(1.0, 1.0)
    @named bound   = PromoterRegion(1.0)
    @named unbound = PromoterRegion(2.0)
    @named model   = RegulatedPromoter(bound, unbound, ligand, 1.0, 1.0)

    @test length(equations(model)) == 2
    @test length(states(model)) == 2
    @test length(ModelingToolkit.get_states(model)) == 0
    @test length(parameters(model)) == 4
    @test length(ModelingToolkit.get_ps(model)) == 2
    @test length(ModelingToolkit.get_systems(model)) == 2
    @test GLG.component_type(ModelingToolkit.get_systems(model)[1]) == GLG.PromoterRegion
    @test GLG.component_type(ModelingToolkit.get_systems(model)[2]) == GLG.PromoterRegion

    @variables t
    @test isequal(t, ModelingToolkit.get_iv(model))

    @parameters k₁ k₀
    @test any(isequal(k₁), parameters(model))
    @test any(isequal(k₀), parameters(model))
    @test ModelingToolkit.get_defaults(model)[k₁] == 1.0
    @test ModelingToolkit.get_defaults(model)[k₀] == 1.0

    @test GLG.component_type(model) == GLG.RegulatedPromoter
end

@testset "Regulated Promoter constructors" begin
    @named ligand  = InputSpecies(1.0, 1.0)
    @named bound   = PromoterRegion(1.0)
    @named unbound = PromoterRegion(2.0)
    @named model   = RegulatedPromoter(bound, unbound, ligand, 1.0, 1.0)

    @test model == RegulatedPromoter(1.0, 2.0, ligand, 1.0, 1.0; name=:model)
    @test model == RegulatedPromoter(bound, unbound, ligand, 1.0, 1.0; name=:model)
end

@testset "Initially random conditions - conservation law for promoters" begin
    @named ligand  = InputSpecies(1.0, 1.0)
    @named bound   = PromoterRegion(1.0)
    @named unbound = PromoterRegion(2.0)
    @named model   = RegulatedPromoter(bound, unbound, ligand, 1.0, 1.0)

    @test sum(last.(collect(randu0(model)))) == 1
end

@testset "InputSpecies definition" begin
    @named model = InputSpecies(1.0, 1.0)

    @test length(equations(model)) == 2
    @test length(states(model)) == 1
    @test length(parameters(model)) == 2

    @variables t species(t)
    @test isequal(t, ModelingToolkit.get_iv(model))
    @test any(isequal(species), states(model))
    @test any(isequal(species), ModelingToolkit.outputs(model))

    @parameters λ α
    @test any(isequal(λ), parameters(model))
    @test any(isequal(α), parameters(model))
    @test ModelingToolkit.get_defaults(model)[λ] == 1.0
    @test ModelingToolkit.get_defaults(model)[α] == 1.0

    @test GLG.component_type(model) == GLG.InputSpecies
end

@testset "Monomer definition" begin
    @named model = Monomer(1.0, 2.0, 3.0)

    @test length(equations(model)) == 3
    @test length(states(model)) == 2
    @test length(parameters(model)) == 3

    @variables t rna(t) monomer(t)
    @test isequal(t, ModelingToolkit.get_iv(model))
    @test any(isequal(rna), states(model))
    @test any(isequal(monomer), states(model))
    @test any(isequal(rna), ModelingToolkit.inputs(model))
    @test any(isequal(monomer), ModelingToolkit.outputs(model))

    @parameters λ α₁ α₂
    @test any(isequal(λ), parameters(model))
    @test any(isequal(α₁), parameters(model))
    @test any(isequal(α₂), parameters(model))
    @test ModelingToolkit.get_defaults(model)[λ] == 1.0
    @test ModelingToolkit.get_defaults(model)[α₁] == 2.0
    @test ModelingToolkit.get_defaults(model)[α₂] == 3.0

    @test GLG.component_type(model) == GLG.Monomer
end

@testset "Dimer definition" begin
    @named model = Dimer(1.0, 2.0, 3.0, 4.0, 5.0)

    @test length(equations(model)) == 5
    @test length(states(model)) == 3
    @test length(parameters(model)) == 5

    @variables t rna(t) monomer(t) dimer(t)
    @test isequal(t, ModelingToolkit.get_iv(model))
    @test any(isequal(rna), states(model))
    @test any(isequal(monomer), states(model))
    @test any(isequal(dimer), states(model))
    @test any(isequal(rna), ModelingToolkit.inputs(model))
    @test any(isequal(dimer), ModelingToolkit.outputs(model))

    @parameters λ α₁ α₂ k₁ k₋₁
    @test any(isequal(λ), parameters(model))
    @test any(isequal(α₁), parameters(model))
    @test any(isequal(α₂), parameters(model))
    @test any(isequal(k₁), parameters(model))
    @test any(isequal(k₋₁), parameters(model))
    @test ModelingToolkit.get_defaults(model)[λ] == 1.0
    @test ModelingToolkit.get_defaults(model)[α₁] == 2.0
    @test ModelingToolkit.get_defaults(model)[α₂] == 3.0
    @test ModelingToolkit.get_defaults(model)[k₁] == 4.0
    @test ModelingToolkit.get_defaults(model)[k₋₁] == 5.0

    @test GLG.component_type(model) == GLG.Dimer
end

@testset "Connections between promoters" begin
    @named input = InputSpecies(1.0, 1.0)
    @named monomer = Monomer(1.0, 1.0, 1.0)
    @named dimer = Dimer(1.0, 1.0, 1.0, 1.0, 1.0)

    @named simplepromoter = PromoterRegion(1.0)
    @named promoter1 = RegulatedPromoter(1.0, 1.0, input, 1.0, 1.0)
    @named promoter2 = RegulatedPromoter(1.0, 1.0, monomer, 1.0, 1.0)
    @named promoter3 = RegulatedPromoter(1.0, 1.0, dimer, 1.0, 1.0)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, monomer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, monomer.rna]),
    ]
    @test rxs == connections(promoter1, promoter2)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, dimer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, dimer.rna]),
    ]
    @test rxs == connections(promoter1, promoter3)
    
    rxs = [
        Reaction(simplepromoter.λ, [simplepromoter.promoter], [simplepromoter.promoter, monomer.rna]),
    ]
    @test rxs == connections(simplepromoter, promoter2)

    rxs = [
        Reaction(simplepromoter.λ, [simplepromoter.promoter], [simplepromoter.promoter, dimer.rna]),
    ]
    @test rxs == connections(simplepromoter, promoter3)

    rxs = [
        Reaction(promoter2.bound.λ, [promoter2.bound.promoter], [promoter2.bound.promoter, monomer.rna]),
        Reaction(promoter2.unbound.λ, [promoter2.unbound.promoter], [promoter2.unbound.promoter, monomer.rna]),
    ]
    @test rxs == connections(promoter2, promoter2)

    rxs = [
        Reaction(promoter1.bound.λ, [promoter1.bound.promoter], [promoter1.bound.promoter, monomer.rna, dimer.rna]),
        Reaction(promoter1.unbound.λ, [promoter1.unbound.promoter], [promoter1.unbound.promoter, monomer.rna, dimer.rna]),
    ]
    @test rxs == connections(promoter1, promoter2, promoter3)
end
