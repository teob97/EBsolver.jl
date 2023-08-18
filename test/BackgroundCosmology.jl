@test BackgroundCosmologySetup(T0_CMB = 0.0).Ω0_γ  ≈ 0.0
@test BackgroundCosmologySetup(N_eff = 0.0).Ω0_nu  ≈ 0.0

# The right sides have been taken from "Cosmology" [Weinberg 2ed]
h_test = 0.67
T_test = 2.7255
@test isapprox(eval_ρ0_crit(h_test*100*(1000/3.086e22)), 1.878e-26 * h_test^2; rtol=1e-4)
@test isapprox(BackgroundCosmologySetup(h = h_test, T0_CMB = T_test).Ω0_γ, 2.47e-5 / h_test^2; rtol=1e-2)

# η'(x) * Hp(x) / c must be ≈ 1
bc = BackgroundCosmology(BackgroundCosmologySetup())
x = range(bc.setup.x_start, bc.setup.x_end, bc.setup.n_splines)
eta_derivative = BSplineKit.diff(bc.η_of_x)
@test isapprox(mean(eta_derivative.(x).*bc.Hp_of_x(x)/c_SI), 1.0, rtol=1e-5)

# Analitically, for a fluid with EoS P=wρ we have Hp'/Hp = d[log(Hp)]/dt = - [(1+3w) / 2] H(t)

# Einstein–de Sitter universe (only matter with w=0) 
bc = BackgroundCosmology(BackgroundCosmologySetup(Ω0_B=0.0, Ω0_CDM=1.0, Ω0_γ=0.0, Ω0_nu=0.0))
@test (bc.dHpdx_of_x(x)./bc.Hp_of_x(x)) ≈ -0.5*bc.H_of_x(x)

# Only radiation w=1/3
bc = BackgroundCosmology(BackgroundCosmologySetup(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=1.0, Ω0_nu=0.0))
@test (bc.dHpdx_of_x(x)./bc.Hp_of_x(x)) ≈ -bc.H_of_x(x)

# Only lambda w=1/3
bc = BackgroundCosmology(BackgroundCosmologySetup(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=0.0, Ω0_nu=0.0))
@test (bc.dHpdx_of_x(x)./bc.Hp_of_x(x)) ≈ bc.H_of_x(x)