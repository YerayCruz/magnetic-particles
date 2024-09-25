module GRAPHIC
using Plots
using Unitful
export graphic
function graphic(x, t, L, σ)
    #x = x |> u"μm"
    #L = L |> u"μm"
    #σ = σ |> u"μm"
    anim = @animate for i = 1:5:length(t)
      
      plot(legend=true, xlims=(-L, L), ylims=(-L, L), xlabel="x", ylabel="y", title="Time: $(round(t[i], digits=2))")
      for b = 1:length(x)
    
          scatter!(x[b]["position"][i, :, 1], x[b]["position"][i, :, 2])
      end
    end
    
    gif(anim, "particle_simulation.gif", fps = 24)
end
end
