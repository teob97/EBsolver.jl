export BackgroundCosmology, eval_ρ0_crit

eval_ρ0_crit(H0) = (3*(H0*(1000/3.086e22))^2)/(8*π*G_SI)
eval_Ω0_γ(T0_CMB_kelvin::Float64, H0::Float64) = 2 * (π^2/30) * (k_b_SI*T0_CMB_kelvin)^4 / (ħ_SI^3 * c_SI^3 * eval_ρ0_crit(H0) * c_SI^2)
eval_Ω0_nu(N_eff::Float64, T0_CMB_kelvin::Float64, H0::Float64) = N_eff * (7.0/8) * (4.0/10)^(4.0/3) * eval_Ω0_γ(T0_CMB_kelvin, H0)
eval_Ω0_Λ(Ω0_k, Ω0_B, Ω0_CDM, Ω0_γ, Ω0_nu)	= 1 - (Ω0_k + Ω0_B + Ω0_CDM + Ω0_γ + Ω0_nu)

Base.@kwdef struct BackgroundCosmology
	
	h       :: Float64 = 0.67
	Ω0_B    :: Float64 = 0.05	 
	Ω0_CDM  :: Float64 = 0.267
	Ω0_k    :: Float64 = 0.0
	T0_CMB  :: Float64 = 2.7255
	N_eff   :: Float64 = 3.046
	
	Ω0_γ    :: Float64 = eval_Ω0_γ(T0_CMB, h*100)
	Ω0_nu   :: Float64 = eval_Ω0_nu(N_eff, T0_CMB, h*100)
	Ω0_Λ    :: Float64 = eval_Ω0_Λ(Ω0_k, Ω0_B, Ω0_CDM, Ω0_γ, Ω0_nu)

	x_start :: Float64 = log(1e-8)
	x_end   :: Float64 = log(1.0)
	n_splines :: Int = 1000

end
