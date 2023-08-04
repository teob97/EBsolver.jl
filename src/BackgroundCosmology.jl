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
	#function DerivedParameters(x::CosmologicalParameters)
	#	Ω0_γ = 2 * (π^2 / 30)
	#end
end
