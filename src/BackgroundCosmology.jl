export BackgroundCosmology 
export eval_ρ0_crit, eval_Ω0_nu, eval_Ω0_γ, eval_Ω0_Λ
export H_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x, η_of_x

eval_ρ0_crit(H0) = (3*H0^2)/(8*π*G_SI)
eval_Ω0_γ(T0_CMB_kelvin::Float64, H0::Float64) = 2 * (π^2/30) * (k_b_SI*T0_CMB_kelvin)^4 / (ħ_SI^3 * c_SI^3 * eval_ρ0_crit(H0) * c_SI^2)
eval_Ω0_nu(N_eff::Float64, T0_CMB_kelvin::Float64, H0::Float64) = N_eff * (7.0/8) * (4.0/10)^(4.0/3) * eval_Ω0_γ(T0_CMB_kelvin, H0)
eval_Ω0_Λ(Ω0_k, Ω0_B, Ω0_CDM, Ω0_γ, Ω0_nu)	= 1 - (Ω0_k + Ω0_B + Ω0_CDM + Ω0_γ + Ω0_nu)

Base.@kwdef struct BackgroundCosmology
	
	h       :: Float64 = 0.67
	H0_SI	:: Float64 = h*100*(1000/3.086e22) #convert Km and Mpc into m
	Ω0_B    :: Float64 = 0.05	 
	Ω0_CDM  :: Float64 = 0.267
	Ω0_k    :: Float64 = 0.0
	T0_CMB  :: Float64 = 2.7255
	N_eff   :: Float64 = 3.046
	
	Ω0_γ    :: Float64 = eval_Ω0_γ(T0_CMB, H0_SI)
	Ω0_nu   :: Float64 = eval_Ω0_nu(N_eff, T0_CMB, H0_SI)
	Ω0_Λ    :: Float64 = eval_Ω0_Λ(Ω0_k, Ω0_B, Ω0_CDM, Ω0_γ, Ω0_nu)

	x_start :: Float64 = log(exp(-12))
	x_end   :: Float64 = log(1.0)
	n_splines :: Int = 1000

end

# We use x=log10(a) where a=a(t) is the scale factor
function H_of_x(x::Float64, BC::BackgroundCosmology)::Float64
	return BC.H0_SI*sqrt((BC.Ω0_B+BC.Ω0_CDM)*exp(-3*x) + (BC.Ω0_γ+BC.Ω0_nu)*exp(-4*x) + BC.Ω0_k*exp(-2*x) + BC.Ω0_Λ)
end

function Hp_of_x(x::Float64, BC::BackgroundCosmology)::Float64
	return exp(x)*H_of_x(x,BC)	
end

function dHpdx_of_x(x::Float64, BC::BackgroundCosmology)::Float64
	return 0.5 * BC.H0_SI^2 * (-(BC.Ω0_B+BC.Ω0_CDM)*exp(-2*x) - 2*(BC.Ω0_γ+BC.Ω0_nu)*exp(-3*x) + 2*BC.Ω0_Λ*exp(x))
end

function ddHpddx_of_x(x::Float64, BC::BackgroundCosmology)::Float64
	return Hp_of_x(x,BC) * 0.5 * BC.H0_SI^2 * (2*(BC.Ω0_B+BC.Ω0_CDM)*exp(-3*x) + 6*(BC.Ω0_γ+BC.Ω0_nu)*exp(-4*x) + 2*BC.Ω0_Λ)
end

function η_of_x(BC::BackgroundCosmology)

	tspan = (BC.x_start, BC.x_end)
	u_0 = c_SI / Hp_of_x(BC.x_start, BC)
	conformal_time(u, p, t) = c_SI / Hp_of_x(t, BC)

	prob = ODE.ODEProblem(conformal_time, u_0, tspan)
	f = ODE.solve(prob)
	
	x = LinRange(BC.x_start, BC.x_end, BC.n_splines)
	return Spline.interpolate(x, f(x).u, Spline.BSplineOrder(3))

end