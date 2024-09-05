module FORCE_CALCULATION

include("./potentials.jl")
using .POTENTIALS

export wall_force, angle_force, create_bond, chain_interaction

function wall_force(x, i, ϵ, σ, b, L, F)

    for ip = 1:x[b].n

      FW_x1        = F_lj_2d([x[b].position[i, ip, 1] - (-L), 0], ϵ, σ)[1]
      FW_x2        = F_lj_2d([x[b].position[i, ip, 1] - (+L), 0], ϵ, σ)[1]
      FW_y1        = F_lj_2d([0, x[b].position[i, ip, 2] - (-L)], ϵ, σ)[2]
      FW_y2        = F_lj_2d([0, x[b].position[i, ip, 2] - (+L)], ϵ, σ)[2]
      FW           = [FW_x1 + FW_x2, FW_y1 + FW_y2]

      F[b].force[i, ip, 1]   -= FW[1]
      F[b].force[i, ip, 2]   -= FW[2]

    end
end


function angle_force(x, i, b, F, ka, angle)
  for ip = 2:2:x[b].n - 1
    r1 = x[b].position[i, ip - 1, :] - x[b].position[i, ip, :]
    r2 = x[b].position[i, ip + 1, :] - x[b].position[i, ip, :]
    angle_force = harmonic_force(r1, r2, angle, ka)
    F[b].force[i, ip - 1, :] += angle_force
  end
end


function create_bond(x, viscocity, v, i, b, F, k, d) #x is an array

   for p = 1:x[b].n - 1

     r = x[b].position[i, p, :] - x[b].position[i, p+1, :]
     elastic_force = elastic_f(r, k, d)
     F[b].force[i, p + 1, :] -= elastic_force
     F[b].force[i, p, : ] += elastic_force
   end
 
end

function chain_interaction(x, b, ϵ, σ, F, numberb, i)
  for nextb = (b + 1):numberb
   for ip = 1:x[b].n
     for jp = 1:x[nextb].n
       r = x[b].position[i, ip, :] - x[nextb].position[i, jp, :]
       F[nextb].force[i, jp, :] += F_lj_2d(r, ϵ, σ)
       F[b].force[i, ip, :] -= F_lj_2d(r, ϵ, σ)
     end
   end
  end
end


end
