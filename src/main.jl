using Plots
using Random
using LinearAlgebra: norm
using Statistics

using ProgressMeter

include("./initial_pos/init_pos.jl")
include("./Forces/forces.jl")

using .BONDS_POSITION, .FORCES

k = 0.1
d = 4.9

bonds = [2, 2, 4]

viscocity = 0.9

kB = 8.617e-5 # Boltzman Constant
ϵ = 1.0
σ = 1.0
wall_dimension = 20
np = 2

T = 300.15 #Temperature of system

m = 1.0 # Mass of particles

L = wall_dimension/2.0

t = range(0, stop=40, length=4001)
struct Velocities
  n::Int
  velocity::Array{Float64, 3}
end

struct Forces
  n::Int
  force::Array{Float64, 3}
end

v = []
F = []

for i = 1:length(bonds)
  n = bonds[i]
  zeros_array = zeros(Float64, length(t), n, 2)

  push!(v, Velocities(n, zeros_array))
end

for i = 1:length(bonds)
  n = bonds[i]
  zeros_array = zeros(Float64, length(t), n, 2)

  push!(F, Forces(n, zeros_array))
end
temperature = zeros(Float64, length(t)) 

x = initial_positions_distribution(σ, L, bonds, 2, t) ## x es un array de n estructuras que contiene el numero de atomos de cada bond y sus posiciones


function wall_force(x, i, ϵ, σ, b)

    for ip = 1:x[b].n

      FW_x1        = F_lj_2d([x[b].position[i, ip, 1] - (-L), 0], ϵ, σ)[1]
      FW_x2        = F_lj_2d([x[b].position[i, ip, 1] - (+L), 0], ϵ, σ)[1]
      FW_y1        = F_lj_2d([0, x[b].position[i, ip, 2] - (-L)], ϵ, σ)[2]
      FW_y2        = F_lj_2d([0, x[b].position[i, ip, 2] - (+L)], ϵ, σ)[2]
      FW           = [FW_x1 + FW_x2, FW_y1 + FW_y2]

      F[b].force[i, ip, 1]   -= FW[1]
      F[b].force[i, ip, 2]   -= FW[2]

      #FW_x1        = F_lj_2d([x[b].position[i, ip + 1, 1] - (-L), 0], ϵ, σ)[1]
      #FW_x2        = F_lj_2d([x[b].position[i, ip + 1, 1] - (+L), 0], ϵ, σ)[1]
      #FW_y1        = F_lj_2d([0, x[b].position[i, ip + 1, 2] - (-L)], ϵ, σ)[2]
      #FW_y2        = F_lj_2d([0, x[b].position[i, ip + 1, 2] - (+L)], ϵ, σ)[2]
      #FW           = [FW_x1 + FW_x2, FW_y1 + FW_y2]

      #F[b].force[i,ip + 1, 1]   += FW[1]
      #F[b].force[i,ip + 1, 2]   += FW[2]

    end
end

#Update force calculation to handle boundaries and particle interactions


function create_bond(x, viscocity, velocity, i, b) #x is an array

   for p = 1:x[b].n - 1

     r = x[b].position[i, p, :] - x[b].position[i, p+1, :]
     elastic_force = elastic_f(r, k, d)
     F[b].force[i, p + 1, :] -= (elastic_force + viscocity * v[b].velocity[i, p + 1, :])
     F[b].force[i, p, : ] += (elastic_force - viscocity * v[b].velocity[i, p, :])
   end
 
end

Δt = step(t)
@showprogress "Computing velocity-verlet" for i in 1:length(t)-1
  for b = 1:length(bonds)

    wall_force(x, i, ϵ, σ, b)
    create_bond(x, viscocity, v, i, b)
    v_mid = v[b].velocity[i, :, :] .+ Δt / 2 / m .* F[b].force[i, :, :]
    x[b].position[i + 1, :, :] = x[b].position[i, :, :] .+ v_mid
    wall_force(x, i + 1, ϵ, σ, b)
    create_bond(x, viscocity, v, i+1, b)
    v[b].velocity[i+1, :, : ] = v_mid .+ Δt / 2 / m .* F[b].force[i+1, :, :]

  end
#    update_forces(x, F, i, ϵ, σ, np, viscocity, v[i, :, :])
#    v_mid = v[i,:,:] .+ Δt / 2 / m .* F[i,:,:]
#    x[i+1,:,:] = x[i,:,:] .+ Δt .* v_mid
#    update_forces(x, F, i+1, ϵ, σ, np, viscocity, v[i, :, :]) # Update forces based on new positions
#    v[i+1,:,:] = v_mid .+ Δt / 2 / m .* F[i+1,:,:]

end


anim = @animate for i = 1:5:length(t)
  
  plot(legend=true, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", title="Time: $(round(t[i], digits=2))", xticks=-10:1:10)
  for b = 1:length(bonds)

    scatter!(x[b].position[i, :, 1], x[b].position[i, :, 2], markersize= σ*10)
  end
end

gif(anim, "particle_simulation.gif", fps = 24)
