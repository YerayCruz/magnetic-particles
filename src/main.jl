using Plots
using ProgressMeter

include("./initial_values/init_structures.jl")
include("./Forces/force_calculation.jl")
include("./Forces/potentials.jl")
using .BONDS_POSITION, .FORCE_CALCULATION

k = 0.5
ka = 0.5
d = 2.7
angles = [180, 100]
bonds = [6, 3]
viscocity = 0.9

ϵ = 1.0
σ = 1.0
wall_dimension = 20

L = wall_dimension/2.0

t = range(0, stop=40, length=4001)

x = initial_positions_distribution(σ, L, bonds, 2, t) ## x es un array de n estructuras que contiene el numero de atomos de cada bond y sus posiciones
v = initial_velocities(bonds, t)
F = initial_forces(bonds, t)

Δt = step(t)
@showprogress "Computing velocity-verlet" for i in 1:length(t)-1
  for b = 1:length(bonds)
    
    wall_force(x, i, ϵ, σ, b, L, F)
    create_bond(x, viscocity, v, i, b, F, k, d)
    angle_force(x, i, b, F, ka, angles[b])
    if length(bonds) > 1
        chain_interaction(x, b, ϵ, σ, F, length(bonds), i)
    end

    x[b].position[i + 1, :, :] = x[b].position[i, :, :] + (Δt / viscocity) * F[b].force[i, :, :]

  end

end


anim = @animate for i = 1:5:length(t)
  
  plot(legend=true, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", title="Time: $(round(t[i], digits=2))", xticks=-10:1:10, yticks=-10:1:10)
  for b = 1:length(bonds)

    scatter!(x[b].position[i, :, 1], x[b].position[i, :, 2], markersize= σ*10)
  end
end

gif(anim, "particle_simulation.gif", fps = 24)
