@component function ConstantSpecies(; name)
    @variables t species(t)
    ReactionSystem(Reaction[], t, [species], []; name=name)
end

random_state(x::ConstantSpecies) = Dict(x.species => rand(0:32))

function set_count(x::ConstantSpecies, v::Int)
    Dict(x.species =>  v)
end

@component function Monomer(; name, λ=1, δ=1, γ=1)
    @parameters λ=λ δ=δ γ=γ
    @reaction_network $name begin
        $λ, rna --> rna + monomer
        $δ, rna --> ∅
        $γ, monomer --> ∅
    end
end

random_state(x::Monomer) = Dict(x.rna => rand(0:32), x.monomer => rand(0:32))

function set_count(x::Monomer, v::Int)
    Dict(x.rna =>  0, x.monomer => v)
end

@component function Dimer(; name, λ=1, δ=1, γ=1, binding=1, unbinding=1)
    @parameters λ=λ δ=δ γ=γ binding=binding unbinding=unbinding
    @reaction_network $name begin
        $λ, rna --> rna + monomer
        $δ, rna --> ∅
        $γ, monomer --> ∅
        $γ, dimer --> ∅
        $(binding), monomer + monomer --> dimer
        $(unbinding), dimer --> monomer + monomer
    end
end

function random_state(x::Dimer)
    Dict(x.rna => rand(0:32), x.monomer => rand(0:16), x.dimer => rand(0:16))
end

function set_count(x::Dimer, v::Int)
    Dict(x.rna =>  0, x.monomer => 0, x.dimer => v)
end
