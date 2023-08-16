using EBsolver
using Plots

let
    bc = BackgroundCosmology()
    x = range(bc.x_start, bc.x_end, bc.n_splines)

    H = Hp_of_x(bc, x)/(1000/3.086e22)/100
    p = plot(
        x, H,
        yscale=:log10, 
        title = "Hp(x) [100 km/s/Mpc]",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/Hp.png")

    eta = η_of_x(bc)
    p = plot(
        x, 
        eta.(x)./3.086e22, 
        yscale=:log10, 
        title = "η(x) [Mpc]",
        xlabel = "x = ln(a)",
        legend = false
        )
    savefig("plots/BackgroundCosmology/eta.png")

    H = Hp_of_x(bc, x)
    p = plot(
        x, 
        eta.(x) .* H ./ c_SI,
        title = "η(x) * Hp(x) / c",
        xlabel = "x = ln(a)",
        legend = false
    )
    savefig("plots/BackgroundCosmology/eta_H.png")

end