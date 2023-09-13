export RecombinationHistory
export Xe_saha_equation, Xe_saha_equation_with_He
export Xe_peebles_equation
export Xe_of_x
export n_e_of_x
export τ_of_x, dτdx_of_x, ddτddx_of_x
export visibility_function_of_x
export sound_horizon

Base.@kwdef struct RecombinationHistory

    cosmo :: BackgroundCosmology = BackgroundCosmology()

    # Helium fraction
    Yp :: Float64 = 0.245

    z_reion :: Float64 = 8.0
    Δz_reion :: Float64 = 0.5
    z_He_reion :: Float64 = 3.5
    Δz_He_reion :: Float64 = 0.5

    # Start and end points for recombination arrays
    x_start :: Float64 = cosmo.x_start
    x_end :: Float64 = cosmo.x_end
    # Numbers of points of Xe,ne array
    npts_rec_arrays :: Int64 = 1500

    # Limit for when to switch between Saha and Peebles
    Xe_saha_limit :: Float64 = 0.99

end

function Xe_saha_equation(RH::RecombinationHistory, x::Float64)
    n_b :: Float64 = RH.cosmo.Ω0_B * eval_ρ0_crit(RH.cosmo.H0_SI) * exp(-3x) / m_H_SI
    T_b :: Float64 = RH.cosmo.T0_CMB * exp(-x)
    k :: Float64 = (ħ_SI^3 * n_b)^(-1) * (m_e_SI*k_b_SI*T_b/2π)^(3/2) * exp(-ϵ0_SI/(k_b_SI*T_b))
    # Avoid subtract huge floats with same magnitude 
    if 4.0/k < 1e-5
        Xe = 1.0 - 1.0/k + 2.0/k^2
    else 
        Xe = 0.5 * (sqrt(k^2+4k) - k)
    end 
    return Xe
end

function Xe_saha_equation_with_He_step(RH::RecombinationHistory, x::Float64, f_e::Float64)
    
    f_e_new = 1.0
    diff = Inf64

    # rhs os Saha's equations
    k1, k2, k3 = 0.0, 0.0, 0.0

    n_b = RH.cosmo.Ω0_B * eval_ρ0_crit(RH.cosmo.H0_SI) * exp(-3x) / m_H_SI
    T_b = RH.cosmo.T0_CMB * exp(-x)
    a = (m_e_SI * k_b_SI * T_b / 2π)^(1.5)

    while diff>1e-8

        k1 = 2.0/(f_e*n_b*ħ_SI^3) * a * exp(-ξ0_SI/(k_b_SI*T_b))
        k2 = 4.0/(f_e*n_b*ħ_SI^3) * a * exp(-ξ1_SI/(k_b_SI*T_b))
        k3 = 1.0/(f_e*n_b*ħ_SI^3) * a * exp(-ϵ0_SI/(k_b_SI*T_b))

        X_H1 = k3 / (1.0 + k3)
        X_He1 = k1 / (1.0 + k1 + k1*k2)
        X_He2 = k1 * k2 / (1.0 + k1 + k1*k2)

        f_e_new = (2X_He2 + X_He1)*RH.Yp/4 + X_H1*(1-RH.Yp)
        diff = abs(f_e_new - f_e)
        f_e = f_e_new

    end

    return f_e/(1-RH.Yp)

end

function Xe_saha_equation_with_He(RH::RecombinationHistory, range::AbstractRange)
    Xe = [Xe_saha_equation_with_He_step(RH, range[begin], 1.0)]
    for i in range
        push!(Xe, Xe_saha_equation_with_He_step(RH, i, last(Xe)))
    end
    return Xe[2:end]
end

function rhs_peebles_equation(RH::RecombinationHistory, Xe, x::Float64)

    H = H_of_x(RH.cosmo, x)

    T_b = RH.cosmo.T0_CMB * exp(-x)

    n_b = RH.cosmo.Ω0_B * eval_ρ0_crit(RH.cosmo.H0_SI) * exp(-3x) / m_H_SI
    n_H = (1 - RH.Yp)*n_b
    n_1s = (1 - Xe)*n_H

    Λ_α = H/8π^2 * (n_1s*ħ_SI^3)^(-1) * (2*ϵ0_SI/c_SI)^2
    Λ_2s_1s = 8.227

    φ_2 = 0.448 * log(ϵ0_SI/(k_b_SI*T_b))
    α2 = 64π/sqrt(27π) * (fine_structure/m_e_SI)^2 * (ħ_SI^2 / c_SI) * sqrt(ϵ0_SI/(k_b_SI*T_b)) * φ_2
    β = α2/ħ_SI^3 * (m_e_SI * k_b_SI * T_b / 2π)^(3/2) * exp(-ϵ0_SI/(k_b_SI*T_b))
    β2 = α2/ħ_SI^3 * (m_e_SI * k_b_SI * T_b / 2π)^(3/2) * exp(-ϵ0_SI/(4*k_b_SI*T_b))

    Cr = (Λ_2s_1s + Λ_α)/(Λ_2s_1s + Λ_α + β2)

    return Cr/H * (β*(1-Xe) - n_H*α2*Xe^2)

end

function Xe_peebles_equation(RH::RecombinationHistory, x_start)

    u_0 = Xe_saha_equation_with_He_step(RH, x_start, 1.0)

	rhs(u, p, t) = rhs_peebles_equation(RH, u, t)

	prob = ODE.ODEProblem(rhs, u_0, (x_start, RH.x_end))
	
	return ODE.solve(prob, ODE.Tsit5(), verbose = false)
    
end

function Xe_no_reion_of_x(RH::RecombinationHistory)

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)
    saha = Xe_saha_equation_with_He(RH, x)

    indx_limit = first(findall(x -> x<RH.Xe_saha_limit, saha))
    x_limit =  x[indx_limit]
    peebles = Xe_peebles_equation(RH,x_limit)(x[indx_limit:end])

    Xe = [saha[begin:indx_limit-1];peebles]
    
    return Spline.interpolate(x, Xe, Spline.BSplineOrder(4))

end

function Xe_reion_of_x(RH::RecombinationHistory)

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)
    saha = Xe_saha_equation_with_He(RH, x)

    indx_limit = first(findall(x -> x<RH.Xe_saha_limit, saha))
    x_limit =  x[indx_limit]
    peebles = Xe_peebles_equation(RH,x_limit)

    # Reionization correction
    f_He = 0.25 * RH.Yp / (1 - RH.Yp)
    y_reion = (1 + RH.z_reion)^(1.5)
    Δy_reion = 1.5 * RH.Δz_reion * sqrt(1+RH.z_reion)
    correction_H = [0.5 * (1+f_He) * (1+tanh((y_reion-exp(-1.5*i))/Δy_reion)) for i in x[indx_limit:end]] 
    correction_He = [0.5 * f_He * (1+tanh((RH.z_He_reion - x_to_redshift(i) )/RH.Δz_He_reion)) for i in x[indx_limit:end]]
    peebles = peebles(x[indx_limit:end]) .+ correction_H .+ correction_He

    Xe = [saha[begin:indx_limit-1];peebles]
    
    return Spline.interpolate(x, Xe, Spline.BSplineOrder(4))

end

function Xe_of_x(RH::RecombinationHistory, reionization::Bool = true)
    if reionization
        return Xe_reion_of_x(RH)
    else
        return Xe_no_reion_of_x(RH)
    end
end

function n_e_of_x(RH::RecombinationHistory, reionization::Bool = true)

    Xe = Xe_of_x(RH, reionization)

    n_b(x) = RH.cosmo.Ω0_B * eval_ρ0_crit(RH.cosmo.H0_SI) * exp(-3x) / m_H_SI
    n_H(x) = (1.0 - RH.Yp) * n_b(x)

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)

    ne = [Xe(i)*n_H(i) for i in x]

    return Spline.interpolate(x, ne, Spline.BSplineOrder(4))

end

function τ_of_x(RH::RecombinationHistory, reionization::Bool = true)

    u_0 = 0.0
    n_e = n_e_of_x(RH, reionization)

	rhs(u, p, t) = - σ_T_SI * c_SI * n_e(t) / H_of_x(RH.cosmo, t)

	prob = ODE.ODEProblem(rhs, u_0, (RH.x_end, RH.x_start))
	
	f = ODE.solve(prob, ODE.Tsit5(), verbose = false)

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)

    return Spline.interpolate(x, f(x), Spline.BSplineOrder(4))

end

function dτdx_of_x(RH::RecombinationHistory, reionization::Bool = true)
    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)
    n_e = n_e_of_x(RH, reionization)
    dτ = [- σ_T_SI * c_SI * n_e(t) / H_of_x(RH.cosmo, t) for t in x]
    return Spline.interpolate(x, dτ, Spline.BSplineOrder(4))
end

function ddτddx_of_x(RH::RecombinationHistory, reionization::Bool = true)
    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)
    dτ = dτdx_of_x(RH, reionization)
    return Spline.interpolate(x, diff(dτ).(x), Spline.BSplineOrder(4))
end

function visibility_function_of_x(RH::RecombinationHistory, reonization::Bool = true)

    τ = τ_of_x(RH, reonization)
    τ_derivative = dτdx_of_x(RH, reonization)

    g(x) = -τ_derivative(x)*exp(-τ(x))

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)

    return Spline.interpolate(x, g.(x), Spline.BSplineOrder(4))

end

function sound_horizon(RH::RecombinationHistory)
    R(x) = (4.0 * RH.cosmo.Ω0_γ) / (3.0 * RH.cosmo.Ω0_B * exp(x))
    c_s(x) = c_SI * sqrt(R(x)/(3.0*(R(x)+1)))

    s_0 = c_s(RH.x_start) / Hp_of_x(RH.cosmo, RH.x_start)

    rhs(u, p, t) = c_s(t) / Hp_of_x(RH.cosmo, t)

	prob = ODE.ODEProblem(rhs, s_0, (RH.x_start, RH.x_end))
	
	f = ODE.solve(prob, ODE.Tsit5(), verbose = false)

    x = range(RH.x_start, RH.x_end, RH.npts_rec_arrays)

    return Spline.interpolate(x, f(x), Spline.BSplineOrder(4))
end
