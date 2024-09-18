module BONDS_POSITION

using Plots
using DynamicQuantities

export initial_positions_distribution, initial_forces, initial_velocities, magnetic_moment

function initial_positions_distribution(sigma, L, particles_per_chain, nd, t)

  chains = []
    for i = 1:length(particles_per_chain)

        n = particles_per_chain[i]
        zeros_array = zeros(Float64, length(t), n, nd)u"nm"

        
        push!(chains, Dict("n" => n, "position" => zeros_array))
    
    end
    
    number_bonds = length(chains)

    bond_spacing = 1.0 * sigma
    particle_spacing = 2.5 * sigma

    bond_delim = number_bonds / 2
    y_o = -bond_spacing * floor(bond_delim)

    if length(chains) * bond_spacing > 2*L - 5*sigma
        println("Error, bonds out of bounds")
        return
    end

    for i = 1 : number_bonds

        if chains[i]["n"] * particle_spacing > 2*L - 5*sigma
            println("Error, particles out of bounds")
            return
        end

        particle_delim = chains[i]["n"] / 2
        x_o = -particle_spacing * floor(particle_delim)

        for j = 1:chains[i]["n"]
            
            chains[i]["position"][1, j, :] = [x_o, y_o]

            x_o += particle_spacing

        end

        y_o += bond_spacing

    end

    return chains

end

function initial_velocities(bonds, t)
  v = []

  for b = 1:length(bonds)
    zeros_array = zeros(Float64, length(t), bonds[b], 2)
    push!(v, Dict("velocity" => zeros_array))
  end
  return v
end

function initial_forces(bonds, t)
  F = []

  for b = 1:length(bonds)
    zeros_array = zeros(Float64, length(t), bonds[b], 2)us"pN"
    push!(F, Dict("force" => zeros_array))
  end

  return F
end

function magnetic_moment(bonds, σ, B)
    χ = 0.4
    m = []
    v = (4 / 3) * π * σ ^ 3 
    μ0 = (4 * π * 10 ^ 5)u"mT*nm/A" 
    for b = 1:length(bonds)
        push!(m, v * χ * B / μ0)
    end
    return m
end

end
