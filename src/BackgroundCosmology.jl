export setup_BackgroundCosmology, BackgroundCosmology
export eval_ρ0_crit
export H_of_x, Hp_of_x, dHpdx_of_x, ddHpddx_of_x
export Ω_B, Ω_CDM, Ω_k, Ω_nu, Ω_γ, Ω_Λ
export η_of_x, t_of_x
export comoving_distance_of_x, angular_distance_of_x, luminosity_distance_of_x

eval_ρ0_crit(H0::Float64) = (3*H0^2)/(8*π*G_SI)
eval_Ω0_γ(T0_CMB_kelvin::Float64, H0::Float64) = 2 * (π^2/30) * (k_b_SI*T0_CMB_kelvin)^4 / (ħ_SI^3 * c_SI^3 * eval_ρ0_crit(H0) * c_SI^2)
eval_Ω0_nu(N_eff::Float64, T0_CMB_kelvin::Float64, H0::Float64) = N_eff * (7.0/8) * (4.0/10)^(4.0/3) * eval_Ω0_γ(T0_CMB_kelvin, H0)
eval_Ω0_Λ(Ω0_k::Float64, Ω0_B::Float64, Ω0_CDM::Float64, Ω0_γ::Float64, Ω0_nu::Float64)	= 1 - (Ω0_k + Ω0_B + Ω0_CDM + Ω0_γ + Ω0_nu)

Base.@kwdef struct setup_BackgroundCosmology
	
	h       :: Float64 = 0.67
	H0_SI	:: Float64 = h*100*(1000/Mpc_to_m(1.0)) #convert Km and Mpc into m
	Ω0_B    :: Float64 = 0.05	 
	Ω0_CDM  :: Float64 = 0.267
	Ω0_k    :: Float64 = 0.0
	T0_CMB  :: Float64 = 2.7255
	N_eff   :: Float64 = 3.046
	
	Ω0_γ    :: Float64 = eval_Ω0_γ(T0_CMB, H0_SI)
	Ω0_nu   :: Float64 = eval_Ω0_nu(N_eff, T0_CMB, H0_SI)
	Ω0_Λ    :: Float64 = eval_Ω0_Λ(Ω0_k, Ω0_B, Ω0_CDM, Ω0_γ, Ω0_nu)

	# We use x=log(a) where a=a(t) is the scale factor
	x_start :: Float64 = log(exp(-12))
	x_end   :: Float64 = log(1.0)
	n_splines :: Int64 = 1000

end

Base.show(io::IO, BC::setup_BackgroundCosmology) = print(
	io, 
	"Info about cosmology class:\n",
	"Ω0_B:\t\t", 	BC.Ω0_B, 	"\n",
	"Ω0_CDM:\t\t", 	BC.Ω0_CDM, 	"\n",
	"Ω0_Λ:\t\t", 	BC.Ω0_Λ, 	"\n",
	"Ω0_K:\t\t", 	BC.Ω0_k, 	"\n",
	"Ω0_nu:\t\t", 	BC.Ω0_nu, 	"\n",
	"Ω0_γ:\t\t", 	BC.Ω0_γ, 	"\n",
	"N_eff:\t\t", 	BC.N_eff, 	"\n",
	"h:\t\t", 		BC.h, 		"\n",
	"T0_CMB:\t\t", 	BC.T0_CMB, 	"\n",
)

H_of_x(BC::setup_BackgroundCosmology, x::Float64)			= BC.H0_SI*sqrt((BC.Ω0_B+BC.Ω0_CDM)*exp(-3x) + (BC.Ω0_γ+BC.Ω0_nu)*exp(-4x) + BC.Ω0_k*exp(-2x) + BC.Ω0_Λ)
H_of_x(BC::setup_BackgroundCosmology, x::Vector{Float64})	= [H_of_x(BC, i) for i in x]
H_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)		= H_of_x(BC, collect(x))

Hp_of_x(BC::setup_BackgroundCosmology, x::Float64)			= exp(x)*H_of_x(BC, x)
Hp_of_x(BC::setup_BackgroundCosmology, x::Vector{Float64})	= [Hp_of_x(BC, i) for i in x]
Hp_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)	= Hp_of_x(BC, collect(x))

dHpdx_of_x(BC::setup_BackgroundCosmology, x::Float64)			= 0.5 * BC.H0_SI^2 * (-(BC.Ω0_B+BC.Ω0_CDM)*exp(-2x) - 2*(BC.Ω0_γ+BC.Ω0_nu)*exp(-3x) + 2*BC.Ω0_Λ*exp(x))
dHpdx_of_x(BC::setup_BackgroundCosmology, x::Vector{Float64})	= [dHpdx_of_x(BC, i) for i in x]
dHpdx_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)		= dHpdx_of_x(BC, collect(x))

ddHpddx_of_x(BC::setup_BackgroundCosmology, x::Float64)			= Hp_of_x(BC, x) * 0.5 * BC.H0_SI^2 * (2*(BC.Ω0_B+BC.Ω0_CDM)*exp(-3x) + 6*(BC.Ω0_γ+BC.Ω0_nu)*exp(-4x) + 2*BC.Ω0_Λ)
ddHpddx_of_x(BC::setup_BackgroundCosmology, x::Vector{Float64})	= [ddHpddx_of_x(BC, i) for i in x]
ddHpddx_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)	= ddHpddx_of_x(BC, collect(x))

Ω_B(BC::setup_BackgroundCosmology, x::Float64)::Float64		= BC.Ω0_B	* BC.H0_SI^2 / (exp(3x) * H_of_x(BC, x)^2)
Ω_CDM(BC::setup_BackgroundCosmology, x::Float64)::Float64	= BC.Ω0_CDM	* BC.H0_SI^2 / (exp(3x) * H_of_x(BC, x)^2)
Ω_γ(BC::setup_BackgroundCosmology, x::Float64)::Float64		= BC.Ω0_γ	* BC.H0_SI^2 / (exp(4x) * H_of_x(BC, x)^2)
Ω_nu(BC::setup_BackgroundCosmology, x::Float64)::Float64	= BC.Ω0_nu	* BC.H0_SI^2 / (exp(4x) * H_of_x(BC, x)^2)
Ω_k(BC::setup_BackgroundCosmology, x::Float64)::Float64		= BC.Ω0_k	* BC.H0_SI^2 / (exp(2x) * H_of_x(BC, x)^2)
Ω_Λ(BC::setup_BackgroundCosmology, x::Float64)::Float64		= BC.Ω0_Λ	* BC.H0_SI^2 /  H_of_x(BC, x)^2

Ω_B(BC::setup_BackgroundCosmology, x::Vector{Float64})		= [Ω_B(BC, i) for i in x]
Ω_CDM(BC::setup_BackgroundCosmology, x::Vector{Float64})	= [Ω_CDM(BC, i) for i in x]
Ω_γ(BC::setup_BackgroundCosmology, x::Vector{Float64})		= [Ω_γ(BC, i) for i in x]
Ω_nu(BC::setup_BackgroundCosmology, x::Vector{Float64})		= [Ω_nu(BC, i) for i in x]
Ω_k(BC::setup_BackgroundCosmology, x::Vector{Float64})		= [Ω_k(BC, i) for i in x]
Ω_Λ(BC::setup_BackgroundCosmology, x::Vector{Float64})		= [Ω_Λ(BC, i) for i in x]

Ω_B(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_B(BC, collect(x))
Ω_CDM(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_CDM(BC, collect(x))	
Ω_γ(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_γ(BC, collect(x))
Ω_nu(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_nu(BC, collect(x))
Ω_k(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_k(BC, collect(x))
Ω_Λ(BC::setup_BackgroundCosmology, x::AbstractRange)		= Ω_Λ(BC, collect(x))

function η_of_x(BC::setup_BackgroundCosmology)::Spline.SplineInterpolation

	u_0 = c_SI / Hp_of_x(BC, BC.x_start)
	conformal_time(u, p, t) = c_SI / Hp_of_x(BC,t)

	prob = ODE.ODEProblem(conformal_time, u_0, (BC.x_start, BC.x_end))
	f = ODE.solve(prob)
	
	x = range(BC.x_start, BC.x_end, BC.n_splines)
	return Spline.interpolate(x, f(x).u, Spline.BSplineOrder(3))

end

function comoving_distance_of_x(BC::setup_BackgroundCosmology, x::Union{Float64, Vector{Float64}})
	η = η_of_x(BC)
	# Rember that we are working with x=log(a) where a is the scale factor
	return η(log(1.0)) .- η.(x)
end
function comoving_distance_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)
	return comoving_distance_of_x(BC, collect(x))
end

function transverse_comoving_distance_of_x(BC::setup_BackgroundCosmology, x::Union{Float64, Vector{Float64}})

	X = comoving_distance_of_x(BC, x)
	k = sqrt(abs(BC.Ω0_k)) * BC.H0_SI * X / c_SI

	BC.Ω0_k < 0 ? r = X .* sin.(k) ./ k :
	BC.Ω0_k > 0 ? r = X .* sinh.(k) ./ k : r = X

	return r

end

function angular_distance_of_x(BC::setup_BackgroundCosmology, x::Union{Float64, Vector{Float64}}) 
	r = transverse_comoving_distance_of_x(BC, x)
	return exp.(x) .* r
end
function angular_distance_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)
	return angular_distance_of_x(BC, collect(x))
end

function luminosity_distance_of_x(BC::setup_BackgroundCosmology, x::Union{Float64, Vector{Float64}})
	r = transverse_comoving_distance_of_x(BC, x)
	return exp.(-x) .* r
end
function luminosity_distance_of_x(BC::setup_BackgroundCosmology, x::AbstractRange)
	return luminosity_distance_of_x(BC, collect(x))
end

# Cosmic time as a function of x=log(a)
function t_of_x(BC::setup_BackgroundCosmology)::Spline.SplineInterpolation

	u_0 = 0.5 / H_of_x(BC, BC.x_start)
	cosmic_time(u, p, t) = 1 / H_of_x(BC,t)

	prob = ODE.ODEProblem(cosmic_time, u_0, (BC.x_start, BC.x_end))
	f = ODE.solve(prob)
	
	x = range(BC.x_start, BC.x_end, BC.n_splines)
	return Spline.interpolate(x, f(x).u, Spline.BSplineOrder(3))

end

struct BackgroundCosmology
	
	setup::setup_BackgroundCosmology
	
	H_of_x::Function
	Hp_of_x::Function
	dHpdx_of_x::Function
	ddHpddx_of_x::Function
	
	Ω_B::Function
	Ω_CDM::Function
	Ω_γ::Function
	Ω_nu::Function
	Ω_k::Function
	Ω_Λ::Function

	comoving_distance_of_x::Function
	angular_distance_of_x::Function
	luminosity_distance_of_x::Function

	η_of_x::Spline.SplineInterpolation
	t_of_x::Spline.SplineInterpolation

	x_equality :: Float64
	x_acceleration :: Float64

	BackgroundCosmology(BC::setup_BackgroundCosmology) = new(
		BC,
		x -> H_of_x(BC, x),
		x -> Hp_of_x(BC, x),
		x -> dHpdx_of_x(BC, x),
		x -> ddHpddx_of_x(BC, x),
		x -> Ω_B(BC, x),
		x -> Ω_CDM(BC, x),
		x -> Ω_γ(BC, x),
		x -> Ω_nu(BC, x),
		x -> Ω_k(BC, x),
		x -> Ω_Λ(BC, x),
		x -> comoving_distance_of_x(BC, x),
		x -> angular_distance_of_x(BC, x),
		x -> luminosity_distance_of_x(BC, x),
		η_of_x(BC),
		t_of_x(BC),
		log((BC.Ω0_γ+BC.Ω0_nu)/(BC.Ω0_CDM+BC.Ω0_B)),
		1.0/3.0 * log((BC.Ω0_CDM+BC.Ω0_B)/BC.Ω0_Λ)
	)

end