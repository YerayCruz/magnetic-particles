using HDF5
using Plots
using Unitful

h5open("msd_simulations.h5", "r") do file
    msd = read(file["273.0 K/1/msd"])
    println(msd)

    t = range(0, stop=5, length=1001)
    test = collect(t)u"s"


    plot(test, msd, label="Averaged MSD", xlabel="Time", ylabel="MSD", title="MSD vs Theory")
    t = plot!(test, msd_theory, label="Theorical MSD", linestyle=:dash)
    savefig(t, "msd.png")
end
