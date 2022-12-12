abstract type Product <: Component end

struct Products <: Component
    reaction_system::ReactionSystem
    Products(x::ReactionSystem) = new(x)
end
function Products(products::Vector{<:Product})
    @variables t
    products = [product.reaction_system for product ∈ products]
    Products(ReactionSystem(Reaction[], t; systems=products, name=:products))
end
function Products(products::Vararg{<:Product})
    @variables t
    products = [product.reaction_system for product ∈ products]
    Products(ReactionSystem(Reaction[], t; systems=products, name=:products))
end

struct ConstantSpecies <: Product
    reaction_system::ReactionSystem
end

struct Monomer <: Product
    reaction_system::ReactionSystem
    Monomer(x::ReactionSystem) = new(x)
end

function Monomer(; name)
    rs = @reaction_network $name begin
        λ, rna --> monomer
        δ, rna --> ∅
        γ, monomer --> ∅
    end λ δ γ
    Monomer(rs)
end
protein(x::Monomer) = x.monomer


struct Dimer <: Product
    reaction_system::ReactionSystem
end
function Dimer(; name)
    rs = @reaction_network $name begin
        λ, rna --> monomer
        δ, rna --> ∅
        γ, monomer --> ∅
        γ, dimer --> ∅
        binding, monomer + monomer --> dimer
        unbinding, dimer --> monomer + monomer
    end λ δ γ binding unbinding
    Dimer(rs)
end
protein(x::Dimer) = x.dimer
