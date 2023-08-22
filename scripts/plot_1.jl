using EBsolver
using Plots
using Logging

let

    @info "Setting up Background Cosmology"

    path = "scripts/plots/BackgroundCosmology/milestone1.png"
    bc = BackgroundCosmology()
    x = range(bc.x_start, bc.x_end, bc.n_splines)

    @info "Making H(x) and η(x) graphs"

    p1 = plot(
        x, 
        Hp_of_x(bc, x)/(1000/Mpc_to_m(1.0))/100,
        yaxis=:log,
        xlim=[bc.x_start, bc.x_end],
        ylim=[10^-1, 10^3],
        legend = false,
        title = "Hp(x) [100 Km/s/Mpc]",
        linewidth=2,
        xlabel="x=log(a)"
    )

    p2 =plot(
        x,
        m_to_Mpc.(η_of_x(bc).(x)),
        yaxis=:log,
        xlim=[bc.x_start, bc.x_end],
        ylim=[10^1, 10^6],
        legend = false,
        title = "η(x) [Mpc]",
        linewidth=2,
        xlabel="x=log(a)"
    )

    p3 = plot(
        x,
        (Hp_of_x(bc, x)).*(η_of_x(bc).(x))/c_SI,
        xlim=[bc.x_start, bc.x_end],
        legend = false,
        title = "η(x)*Hp(x)/c",
        linewidth=2,
        xlabel="x=log(a)"
    )

    x1 = range(-20, 5, bc.n_splines)
    p4 = plot(x1, Ω_B(bc, x1)+Ω_CDM(bc, x1), label = "Matter", style=:dash, legend=:left, linewidth=2, xlabel="x=log(a)")
    plot!(x1, Ω_γ(bc, x1)+Ω_nu(bc, x1), label = "Radiation", style=:dashdot, linewidth=2,)
    plot!(x1, Ω_Λ(bc, x1), label = "Lambda", linewidth=2,)

    @info "Start mcmc simulation"

    # h, omega_cdm, omega_k
    data = mcmc_fit_to_supernova_data("scripts\\data\\sn.txt")
    x = data[3]

    mask_1 =  x[1,:] .- data[2] .< 3.53
    mask_2 =  x[1,:] .- data[2] .< 8.02
    mask_3 =  x[1,:] .- data[2] .< 14.16

    @info "Making mcmc graphs"

    p6 = histogram(x[2,mask_3], xlabel="h", normalize=true, label = false)
    p6 = vline!([0.67], linewidth=3, color=:black, style=:dash, label = "Planck best fit")

    omega_r = [eval_Ω0_γ(2.7255, i) for i in x[2,:]]
    lambda = 1 .- (0.05 .+ x[3,:] .+ x[4,:] .+ omega_r)

    p5 = scatter(x[3,mask_2], lambda[mask_2], label = "2σ constrains", legend=:right)
    scatter!(x[3,mask_1], lambda[mask_1], label = "1σ constrains")
    plot!([1.0, 0.0], [0.0,1.0], linewidth = 2, label = "Flat Universe", style=:dash, color=:black)
    xaxis!("Ω_CDM")
    yaxis!("Ω_Λ")

    plot(p1, p2, p3, p4, p5, p6, size = (1000, 600))
    savefig(path)

    @info "Image saved at "*path

end