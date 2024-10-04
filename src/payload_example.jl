include("./brownian_simulation.jl")
include("./difussion.jl")
using .BROWNIAN_SIMULATION
using .DIFUSSION
using Statistics
using Plots

payload = Dict("k" =>  500.0, ##pN/nm
"ka"  => 1000.0, #pN/°
"ϵ" =>  1000.0, #pN*nm
"σ"  => 10.0, #μm
"d" =>  20.7, #μm
"angles" =>  [180, 100], #°
"bonds" =>  [6, 3],
"B" => [0.0, 1.0], #mT
"viscocity" =>  9.0, #pN*s/nm
"wall_dimension" =>  200.0, #μm
"T" => 273.0 #K
)

#brownian_motion(payload)
for T = 273.0:1.0:300.0
    msd_group = []
    for y = 1:800

        payload["T"] = T
        msd, test, msd_theory = difussion(payload, y)

        push!(msd_group, msd)
        if length(msd_group) == 800
            mean_vector = mean(hcat(msd_group...), dims=2)[:]

            plot(test, mean_vector, label="Averaged MSD", xlabel="Time", ylabel="MSD", title="MSD vs Theory $(T)")
            g = plot!(test, 2 .* msd_theory, label="Theorical MSD", linestyle=:dash)
            savefig(g, "~/Documents/Projects/magnetic-particles/mean_msd/mean_msd_$(T).png")
        end

    end
end
