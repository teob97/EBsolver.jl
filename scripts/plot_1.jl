using EBsolver
using Plots

let
    bc = BackgroundCosmology()
    x = LinRange(bc.x_start, bc.x_end, bc.n_splines)

    H = [Hp_of_x(i, bc)/(1000/3.086e22)/100 for i in x]
    p = plot(
        x, H,
        yscale=:log10, 
        title = "Hp(x) [100 km/s/Mpc]",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/Hp.png")

    eta = η_of_x(bc)(x)
    p = plot(
        eta.t, 
        eta.u/3.086e22, 
        yscale=:log10, 
        title = "η(x) [Mpc]",
        xlabel = "x = ln(a)",
        legend = false
        )
    savefig("plots/BackgroundCosmology/eta.png")

    H = [Hp_of_x(i, bc) for i in x]
    p = plot(
        eta.t, 
        (eta.u) .* H / c_SI,
        title = "η(x) * Hp(x) / c",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/eta_H.png")

end