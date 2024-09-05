module BONDS_POSITION

using Plots

export initial_positions_distribution, initial_forces, initial_velocities

struct Chain
    n::Int
    position::Array{Float64, 3}
end


struct Velocities
  velocity::Array{Float64}
end


struct Forces
  force::Array{Float64}
end

#At the end we have an array called chains with structures

function initial_positions_distribution(sigma, L, particles_per_chain, nd, t)

  chains = []
    for i = 1:length(particles_per_chain)

        n = particles_per_chain[i]
        zeros_array = zeros(Float64, length(t), n, nd)
    
        push!(chains, Chain(n, zeros_array))
    
    end
    
    number_bonds = length(chains)

    bond_spacing = 2.5 * sigma
    particle_spacing = 2.5 * sigma

    bond_delim = number_bonds / 2
    y_o = -bond_spacing * floor(bond_delim)

    if length(chains) * bond_spacing > 2*L - 5*sigma
        println("Error, bonds out of bounds")
        return
    end

    for i = 1 : number_bonds

        if chains[i].n * particle_spacing > 2*L - 5*sigma
            println("Error, particles out of bounds")
            return
        end

        particle_delim = chains[i].n / 2
        x_o = -particle_spacing * floor(particle_delim)

        for j = 1:chains[i].n
            
            chains[i].position[1, j, :] = [x_o, y_o]

            x_o += particle_spacing

        end

        y_o += bond_spacing

    end

    plot(xlims = (-L, L), ylims = (-L, L))
    for i = 1:number_bonds
    
        scatter!([chains[i].position[1, :, 1]], [chains[i].position[1, :, 2]])

    end

    savefig("positions.png")

    return chains

end

function initial_velocities(bonds, t)
  v = []

  for b = 1:length(bonds)
    zeros_array = zeros(Float64, length(t), bonds[b], 2)
    push!(v, Velocities(zeros_array))
  end
  return v
end
function initial_forces(bonds, t)
  F = []

  for b = 1:length(bonds)
    zeros_array = zeros(Float64, length(t), bonds[b], 2)
    push!(F, Forces(zeros_array))
  end

  return F
end


end
