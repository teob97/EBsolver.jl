using EBsolver
using Plots

let
    bc = BackgroundCosmology()
    x = LinRange(bc.x_start, bc.x_end, bc.n_splines)

    H = bc.Hp_of_x.(x)/(1000/3.086e22)/100
    p = plot(
        x, H,
        yscale=:log10, 
        title = "Hp(x) [100 km/s/Mpc]",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/Hp.png")

    eta = bc.η_of_x.(x)
    p = plot(
        x, 
        eta/3.086e22, 
        yscale=:log10, 
        title = "η(x) [Mpc]",
        xlabel = "x = ln(a)",
        legend = false
        )
    savefig("plots/BackgroundCosmology/eta.png")

    H = bc.Hp_of_x.(x)
    p = plot(
        x, 
        (eta) .* H / c_SI,
        title = "η(x) * Hp(x) / c",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/eta_H.png")

end