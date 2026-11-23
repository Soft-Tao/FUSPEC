"""A two-reactant, two-product reaction description."""
struct Reaction
    name::Symbol
    reactant1::Particle
    reactant2::Particle
    product1::Particle
    product2::Particle
    link_to_product1::Int8
    link_to_product2::Int8

    function Reaction(name::Symbol, reactant1::Particle, reactant2::Particle,
                      product1::Particle, product2::Particle;
                      link_to_product1::Integer=0,
                      link_to_product2::Integer=0)
        link_to_product1 in 0:2 || throw(ArgumentError(
            "link_to_product1 must be 0, 1, or 2",
        ))
        link_to_product2 in 0:2 || throw(ArgumentError(
            "link_to_product2 must be 0, 1, or 2",
        ))
        new(name, reactant1, reactant2, product1, product2,
            Int8(link_to_product1), Int8(link_to_product2))
    end
end
