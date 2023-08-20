@test BackgroundCosmology(T0_CMB = 0.0).Ω0_γ  ≈ 0.0
@test BackgroundCosmology(N_eff = 0.0).Ω0_nu  ≈ 0.0

# The right sides have been taken from "Cosmology" [Weinberg 2ed]
h_test = 0.67
T_test = 2.7255
@test isapprox(eval_ρ0_crit(h_test*100*(1000/3.086e22)), 1.878e-26 * h_test^2; rtol=1e-4)
@test isapprox(BackgroundCosmology(h = h_test, T0_CMB = T_test).Ω0_γ, 2.47e-5 / h_test^2; rtol=1e-2)

# Test cosmological parameters


# η'(x) * Hp(x) / c must be ≈ 1
bc = BackgroundCosmology()
x = range(bc.x_start, bc.x_end, bc.n_splines)
eta_derivative = BSplineKit.diff(η_of_x(bc))
@test isapprox(mean(eta_derivative.(x).*Hp_of_x(bc, x)/c_SI), 1.0, rtol=1e-5)

# Analitically, for a fluid with EoS P=wρ we have Hp'/Hp = d[log(Hp)]/dt = - [(1+3w) / 2] H(t)
# Einstein–de Sitter universe (only matter with w=0) 
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=1.0, Ω0_γ=0.0, Ω0_nu=0.0)
@test (dHpdx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ -0.5*H_of_x(bc, x)
# Only radiation w=1/3
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=1.0, Ω0_nu=0.0)
@test (dHpdx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ -H_of_x(bc, x)
# Only lambda w=1/3
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=0.0, Ω0_nu=0.0)
@test (dHpdx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ H_of_x(bc, x)

# Analitically, for a fluid with EoS P=wρ we have Hp''/Hp = Hp'/Hp d[log(Hp')]/dt = [(2+9w+9w²) / 2] H(t)
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=1.0, Ω0_γ=0.0, Ω0_nu=0.0)
@test (ddHpddx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ H_of_x(bc, x).^2
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=1.0, Ω0_nu=0.0)
@test (ddHpddx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ 3*H_of_x(bc, x).^2
bc = BackgroundCosmology(Ω0_B=0.0, Ω0_CDM=0.0, Ω0_γ=0.0, Ω0_nu=0.0)
@test (ddHpddx_of_x(bc, x)./Hp_of_x(bc, x)) ≈ H_of_x(bc, x).^2
