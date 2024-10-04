module BROWNIAN_SIMULATION

using ProgressMeter
using Distributions
using Unitful
include("./initial_values/init_structures.jl")
include("./Forces/force_calculation.jl")
include("./Forces/potentials.jl")
include("./graphics.jl")
using .BONDS_POSITION, .FORCE_CALCULATION, .GRAPHIC

export brownian_motion

#@register_unit mT 1e-12u"T"
#@register_unit pN 1e-12u"N"

## Forma del diccionario
#contantes = {
#k = consante de resorte,
#ka - constante angular
#ϵ = consante de interaccion WCA
#σ = radio de las particulas
#d = distancia de equilibrio resorte
#angles = Arreglo de angulos por cadena (proximamente sera por triplete de particulas de cada cadena)
#bonds = Arreglo de numero de particulas por cadena
#B = Arreglo del campo magnetico externo
#viscodity = viscocidad del sistema
#wall dimension - dimensiones de la caja
#T temperature in K
#}


function brownian_motion(payload, show_graphic=true)
    k = payload["k"]u"pN/nm"
    ka = payload["ka"]u"pN"
    ϵ = payload["ϵ"]u"pN*nm"
    σ = payload["σ"]u"μm" |> u"nm"
    d = payload["d"]u"μm" |> u"nm"
    angles = payload["angles"]
    bonds = payload["bonds"]
    B = payload["B"]u"mT"
    viscocity = payload["viscocity"]u"pN*s/nm"
    wall_dimension = payload["wall_dimension"]u"μm" |> u"nm"
    T = payload["T"]u"K"
    kb = 0.013806u"pN*nm/K"

    L = wall_dimension/2.0
    t = range(0, stop=5, length=1001)
    
    x = initial_positions_distribution(σ, L, bonds, 2, t) ## x es un array de n diccionarios que contiene el numero de atomos de cada bond y sus posiciones
    F = initial_forces(bonds, t)
    μ = 1 / viscocity

    m = magnetic_moment(bonds, σ, B)
    Δt = step(t)u"s"

    @showprogress "Computing velocity-verlet" for i in 1:length(t)-1
      for b = 1:length(bonds)
        
        wall_force(x, i, ϵ, σ, b, L, F)
        create_bond(x, viscocity, i, b, F, k, d)
        angle_force(x, i, b, F, ka, angles[b])
        if length(bonds) > 1
            magnetic_interaction(x, b, m, F, length(bonds), i)
            chain_interaction(x, b, ϵ, σ, F, length(bonds), i)
        end

        x[b]["position"][i + 1, :, :] = x[b]["position"][i, :, :] .+ (Δt * μ ) .* F[b]["force"][i, :, :] .+ sqrt(Δt * 2 * kb * T * μ) * rand(Normal())
      end
    end
    
    if show_graphic == true
        graphic(x, t, L, σ)
    end
end

end
