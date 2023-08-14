export BackgroundCosmology 
export eval_ρ0_crit, eval_Ω0_nu, eval_Ω0_γ, eval_Ω0_Λ

eval_ρ0_crit(H0::Float64) = (3*H0^2)/(8*π*G_SI)
eval_Ω0_γ(T0_CMB_kelvin::Float64, H0::Float64) = 2 * (π^2/30) * (k_b_SI*T0_CMB_kelvin)^4 / (ħ_SI^3 * c_SI^3 * eval_ρ0_crit(H0) * c_SI^2)
eval_Ω0_nu(N_eff::Float64, T0_CMB_kelvin::Float64, H0::Float64) = N_eff * (7.0/8) * (4.0/10)^(4.0/3) * eval_Ω0_γ(T0_CMB_kelvin, H0)
eval_Ω0_Λ(Ω0_k::Float64, Ω0_B::Float64, Ω0_CDM::Float64, Ω0_γ::Float64, Ω0_nu::Float64)	= 1 - (Ω0_k + Ω0_B + Ω0_CDM + Ω0_γ + Ω0_nu)

function eval_η_of_x(
	Hp_of_x::Function,
	x_start::Float64, 
	x_end::Float64, 
	n_splines::Int64,
	)::Spline.SplineInterpolation

	u_0 = c_SI / Hp_of_x(x_start)
	conformal_time(u, p, t) = c_SI / Hp_of_x(t)

	prob = ODE.ODEProblem(conformal_time, u_0, (x_start, x_end))
	f = ODE.solve(prob)
	
	x = LinRange(x_start, x_end, n_splines)
	return Spline.interpolate(x, f(x).u, Spline.BSplineOrder(3))

end

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
	n_splines :: Int64 = 1000

	# We use x=ln(a) where a=a(t) is the scale factor
	H_of_x		:: Function	= x -> H0_SI*sqrt((Ω0_B+Ω0_CDM)*exp(-3x) + (Ω0_γ+Ω0_nu)*exp(-4x) + Ω0_k*exp(-2x) + Ω0_Λ)
	Hp_of_x		:: Function	= x -> exp(x)*H_of_x(x)	
	dHpdx_of_x	:: Function	= x -> 0.5 * H0_SI^2 * (-(Ω0_B+Ω0_CDM)*exp(-2x) - 2*(Ω0_γ+Ω0_nu)*exp(-3x) + 2*Ω0_Λ*exp(x))
	ddHpddx_of_x:: Function	= x -> Hp_of_x(x) * 0.5 * H0_SI^2 * (2*(Ω0_B+Ω0_CDM)*exp(-3x) + 6*(Ω0_γ+Ω0_nu)*exp(-4x) + 2Ω0_Λ)

	Ω_B		:: Function = x -> Ω0_B * H0_SI^2 / (exp(3x) * H_of_x(x))
	Ω_CDM	:: Function = x -> Ω0_CDM * H0_SI^2 / (exp(3x) * H_of_x(x))
	Ω_γ		:: Function = x -> Ω0_γ * H0_SI^2 / (exp(4x) * H_of_x(x))
	Ω_nu	:: Function = x -> Ω0_nu * H0_SI^2 / (exp(4x) * H_of_x(x))
	Ω_k		:: Function = x -> Ω0_k * H0_SI^2 / (exp(2x) * H_of_x(x))
	Ω_Λ		:: Function = x -> Ω0_Λ * H0_SI^2 /  H_of_x(x)

	η_of_x :: Spline.SplineInterpolation = eval_η_of_x(Hp_of_x, x_start, x_end, n_splines)

end
