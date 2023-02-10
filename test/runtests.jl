using GeneticLogicGraph
using Test
using Random

Random.seed!(7)

@testset verbose=true "GeneticLogicGraph.jl" begin
    @testset "Graph utilities" begin
        include("utils.jl")
    end

    @testset "Model definitions" begin
        include("models.jl")
    end
end
