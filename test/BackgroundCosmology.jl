@test BackgroundCosmology(T0_CMB = 0.0).Ω0_γ  ≈ 0.0
@test BackgroundCosmology(N_eff = 0.0).Ω0_nu  ≈ 0.0

# The right sides have been taken from "Cosmology" [Weinberg 2ed]
h_test = 0.67
T_test = 2.7255
@test isapprox(eval_ρ0_crit(h_test*100*(1000/3.086e22)), 1.878e-26 * h_test^2; rtol=1e-4)
@test isapprox(BackgroundCosmology(h = h_test, T0_CMB = T_test).Ω0_γ, 2.47e-5 / h_test^2; rtol=1e-2)

@testset begin
    bc = BackgroundCosmology()
    x = LinRange(bc.x_start, bc.x_end, bc.n_splines)
    H = [Hp_of_x(i, bc)/(1000/3.086e22)/100 for i in x]
    p = plot(x, H, yscale=:log10, title = "Hp(x) [100 km/s/Mpc]")
    savefig("../plots/prova_Hp.png")
end

@testset begin
    bc = BackgroundCosmology()
    x = LinRange(bc.x_start, bc.x_end, bc.n_splines)
    f = η_of_x(bc)(x)
    p = plot(f.t, f.u/3.086e22, yscale=:log10, title = "η(x) [Mpc]")
    savefig("../plots/prova_eta.png")
end