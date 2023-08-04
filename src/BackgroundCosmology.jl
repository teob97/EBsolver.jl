Base.@kwdef struct CosmologicalParameters
	h		:: Float64 = 0.67
	Ω0_B	:: Float64 = 0.05	 
	Ω0_CDM	:: Float64 = 0.267
	Ω0_k	:: Float64 = 0.0
	T0_CMB 	:: Float64 = 2.7255
	N_eff 	:: Float64 = 3.046
end

struct DerivedParameters
	Ω0_γ	:: Float64
	Ω0_nu	:: Float64
	Ω0_Λ	:: Float64
	# Constuctor
	function DerivedParameters(x::CosmologicalParameters)
		H0 		= 100 * x.h
		Ω0_γ 	= (16 * π^3 * G_SI * (k_b_SI * x.T0_CMB)^4) / (90 * ħ_SI^3 * c_SI^5 * H0^2)
		Ω0_nu 	= x.N_eff * (7.0/8) * (4.0/10)^(4.0/3) * Ω0_γ
		Ω0_Λ	= 1 - (x.Ω0_k + x.Ω0_B + x.Ω0_CDM + Ω0_γ + Ω0_nu)
		new(Ω0_γ, Ω0_nu, Ω0_Λ)
	end
end
