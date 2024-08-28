using Plots
using Random
using LinearAlgebra: norm
using Statistics

using ProgressMeter

k = 0.1
d = 5

viscocity = 0.7

kB = 8.617e-5 # Boltzman Constant
ϵ = 1.0
σ = 1.0
wall_dimension = 20
np = 2

T = 300.15 #Temperature of system

m = 1.0 # Mass of particles

L = wall_dimension/2.0

t = range(0, stop=40, length=4001)
x = zeros(Float64, length(t), np, 2)
v = zeros(Float64, length(t), np, 2)
F = zeros(Float64, length(t), np, 2)
temperature = zeros(Float64, length(t)) 


# Generate unique initial conditions
function initial_positions_distribution()
    p_dist = 1.5*σ
    x = zeros(Float64, np, 2)
    first_x = -L + p_dist
    first_y = -L + p_dist
    x[1,1] = first_x
    x[1,2] = first_y
    for i=2:np
        if first_x < L - p_dist
            first_x += p_dist
            x[i,1] = first_x
            x[i,2] = first_y
        else
            first_x = -L + p_dist
            first_y += p_dist
            if first_y > L - p_dist
                println("No more space in system")
                return x[1:i-1, :] # Return already filled positions if out of space
            end
            x[i,1] = first_x
            x[i,2] = first_y
        end
    end
    return x
end

x[1,:,:] = initial_positions_distribution()

function elastic_f(r, k, d)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm
    magnitude = -k*(r_norm-d)
    return force_direction * magnitude
end

function F_lj_2d(r, ϵ, σ)
    r_norm = sqrt(sum(r .^ 2))
    force_direction = r / r_norm
    magnitude = 4*ϵ * (12 * (σ / r_norm) ^ 12 - 6 * (σ / r_norm) ^ 6) / r_norm
    return force_direction * magnitude
end

# Update force calculation to handle boundaries and particle interactions
function update_forces(x, F, i, ϵ, σ, np, viscocity, velocity)
    j = 2

    for ip = 1:np - 1
        
        FW_x1        = F_lj_2d([x[i, ip, 1].-(-L), 0], ϵ, σ)[1]
        FW_x2        = F_lj_2d([x[i, ip, 1].-(+L), 0], ϵ, σ)[1]
        FW_y1        = F_lj_2d([0, x[i, ip, 2].-(-L)], ϵ, σ)[2]
        FW_y2        = F_lj_2d([0, x[i, ip, 2].-(+L)], ϵ, σ)[2]
        FW           = [FW_x1 + FW_x2, FW_y1 + FW_y2]
        F[i,ip, 1]   += FW[1]
        F[i,ip, 2]   += FW[2]

        FW_x1        = F_lj_2d([x[i, j, 1].-(-L), 0], ϵ, σ)[1]
        FW_x2        = F_lj_2d([x[i, j, 1].-(+L), 0], ϵ, σ)[1]
        FW_y1        = F_lj_2d([0, x[i, j, 2].-(-L)], ϵ, σ)[2]
        FW_y2        = F_lj_2d([0, x[i, j, 2].-(+L)], ϵ, σ)[2]
        FW           = [FW_x1 + FW_x2, FW_y1 + FW_y2]
        F[i,j, 1]   += FW[1]
        F[i,j, 2]   += FW[2]

        for jp=j:np
            r = x[i, ip, :] - x[i, jp, :]
            force = elastic_f(r, k, d)
            F[i, ip, :] += force - viscocity*velocity[ip, :]
            F[i, jp, :] -= force + viscocity*velocity[jp, :]
        end
        j += 1
    end
end

Δt = step(t)
@showprogress "Computing velocity-verlet" for i in 1:length(t)-1
    
    update_forces(x, F, i, ϵ, σ, np, viscocity, v[i, :, :])
    v_mid = v[i,:,:] .+ Δt / 2 / m .* F[i,:,:]
    x[i+1,:,:] = x[i,:,:] .+ Δt .* v_mid
    update_forces(x, F, i+1, ϵ, σ, np, viscocity, v[i, :, :]) # Update forces based on new positions
    v[i+1,:,:] = v_mid .+ Δt / 2 / m .* F[i+1,:,:]

end


anim = @animate for i = 1:5:length(t)
    plot(legend=false, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", title="Time: $(round(t[i], digits=2))", xticks=-10:1:10)
    for j = 1:np
        scatter!([x[i,j,1]], [x[i,j,2]], markersize= σ*10)
    end
end

gif(anim, "particle_simulation.gif", fps = 24)