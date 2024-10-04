module DIFUSSION
using Unitful 
using ProgressMeter
using Distributions
using Plots
using HDF5

export difussion

function difussion(payload, y)
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

    μ = 1 / viscocity

    L = wall_dimension/2.0
    t = range(0, stop=5, length=1001)u"s"

    n = 1
    x = zeros(Float64, length(t), n, 2)u"nm"
    x[1, 1, :] = [0., 0.]u"nm"
    t = range(0, stop=5, length=1001)
    test = collect(t)u"s"

    Δt = step(t)u"s"

    D = kb * T / viscocity 
    msd_theory = 2 * D * test  # Theoretical MSD for 2D Brownian motion

    @showprogress "Computing velocity-verlet" for i in 1:length(t) - 1
        x[i + 1, 1, :] = x[i, 1, :] .+ sqrt(Δt * 2 * kb * T * μ) * rand(Normal())
    end
    
    #anim = @animate for i = 1:length(t)
    #    plot(legend=true, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", title="Time: $(round(t[i], digits=2))")
    #    scatter!([x[i, 1]], [x[i, 2]])
    #end
    #p = plot(x[:, 1, 1], x[:, 1, 2], legend=true, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", show=true)#, title="Time: $(round(t[i], digits=2))")
    #gif(anim, "difussion.gif", fps=24)
    #savefig(p, "difussion.png")

    N = length(t)
    n = size(x[1, :, :], 1)
    msd = zeros(N)u"nm^2"

    for i = 1:N

        diff = x[i, 1, :] - x[1, 1, :]
        msd[i] = (1 / n) * sum(diff.^2)

    end

    return msd, test, msd_theory
#    x = x ./ 1000.0
#    L = L ./ 1000.0
#
#    plot(test, msd, label="Averaged MSD", xlabel="Time", ylabel="MSD", title="MSD vs Theory")
#    t = plot!(test, msd_theory, label="Theorical MSD", linestyle=:dash)
#    savefig(t, "msd.png")
end
end

