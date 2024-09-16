module MAIN_SIMULATION

using ProgressMeter
include("./initial_values/init_structures.jl")
include("./Forces/force_calculation.jl")
include("./Forces/potentials.jl")
include("./graphics.jl")
using .BONDS_POSITION, .FORCE_CALCULATION, .GRAPHIC

export brownian_motion

## Forma del diccionario
#contantes = {
#k = consante de resorte,
#ka - constante angular
#ϵ = consante de interaccion WCA
#σ = radio de las particulas
#d = distancia de equilibrio resorte
#angles = Arreglo de angulos por cadena (proximamente sera por triplete de particulas de cada cadena)
#bonds = Arreglo de numero de particulas por cadena
#m = Arreglo de los momentos magneticos de cada cadena
#viscodity = viscocidad del sistema
#wall dimension - dimensiones de la caja
#}


function brownian_motion(payload, show_graphic=true)
    k = payload["k"]
    ka = payload["ka"]
    ϵ = payload["ϵ"]
    σ = payload["σ"]
    d = payload["d"]
    angles = payload["angles"]
    bonds = payload["bonds"]
    m = payload["m"]
    viscocity = payload["viscocity"]
    wall_dimension = payload["wall_dimension"]

    L = wall_dimension/2.0
    t = range(0, stop=5, length=1001)
    
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
            magnetic_interaction(x, b, m, F, length(bonds), i)
            chain_interaction(x, b, ϵ, σ, F, length(bonds), i)
        end
    
        x[b].position[i + 1, :, :] = x[b].position[i, :, :] + (Δt / viscocity) .* F[b].force[i, :, :]
    
      end
    
    end
    
    if show_graphic == true
        graphic(x, t, L, σ)
    end
end

end
