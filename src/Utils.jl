export G_SI, c_SI, k_b_SI, k_e_SI, ħ_SI, eV_SI, m_e_SI, m_H_SI, Msun_SI, σ_T_SI, λ_2s1s_SI, pc_SI, fine_structure
export ϵ0_eV, ϵ0_SI, ξ0_eV, ξ1_eV, ξ0_SI, ξ1_SI
export m_to_pc, m_to_Mpc, m_to_Gpc
export pc_to_m, Mpc_to_m, Gpc_to_m
export s_to_Gyr
export x_to_redshift, redshift_to_x
export mcmc_fit_to_supernova_data

# SI constants
const G_SI      = 6.67430e-11
const c_SI      = 2.99792458e8 
const k_b_SI    = 1.38064852e-23
const k_e_SI    = 8.9875517923e9
const ħ_SI      = 1.054571817e-34
const eV_SI     = 1.60217653e-19
const m_e_SI    = 9.10938356e-31
const m_H_SI    = 1.6735575e-27
const Msun_SI   = 1.98847e30
const σ_T_SI    = 6.6524587158e-29
const λ_2s1s_SI = 8.227
const ϵ0_eV     = 13.605693122994
const ϵ0_SI     = ϵ0_eV * eV_SI
const ξ0_eV     = 24.587387
const ξ0_SI     = ξ0_eV * eV_SI
const ξ1_eV     = 4.0 * ξ0_eV
const ξ1_SI     = ξ1_eV * eV_SI
const pc_SI     = 3.08567758e16
const fine_structure = 1.0/137.0359992

# Conversions
m_to_pc(x)  = x * 3.08567758130573e-16
m_to_Mpc(x) = x * 3.08567758130573e-22
m_to_Gpc(x) = x * 3.08567758130573e-26

pc_to_m(x)  = x * 3.08567758130573e16
Mpc_to_m(x) = x * 3.08567758130573e22
Gpc_to_m(x) = x * 3.08567758130573e26

s_to_Gyr(x) = x / (1e9*365*24*60*60)

x_to_redshift(x) = expm1(-x)
redshift_to_x(x) = - log1p(x)

# Simple Monte Carlo Markov Chain with Metropolis sampling to fit supernova redshift vs luminosity distance 

function likelihood_chi2(data::Vector, parameters::Vector{Float64})

    z = [i[1] for i in data]
    x = redshift_to_x.(z)
    d_l_obs = [i[2] for i in data]
    sigma = [i[3] for i in data]

    Ω0_B = 0.05    
    N_eff = 0.00

    setup = BackgroundCosmology(
        h = parameters[1],
        Ω0_B = Ω0_B,
        Ω0_CDM = parameters[2] - Ω0_B,
        Ω0_k = parameters[3],
        N_eff = N_eff
    )

    d_l = luminosity_distance_of_x(setup, x)
    d_l = m_to_Gpc.(d_l)

    return  sum((d_l - d_l_obs).^2 ./ sigma.^2)
end

function mcmc_fit_to_supernova_data(path_to_sn_data::String; likelihood::Function = likelihood_chi2, max_step::Int = 20000)

    # Load data
    file_stream = open(path_to_sn_data, "r")
    readline(file_stream)
    sn_data = []
    while !eof(file_stream)
        x = parse.(Float64, split(readline(file_stream)))
        push!(sn_data, x)
    end
    close(file_stream)

    # Set prior for h, Ω_CDM, Ω_k
    prior_upper = [1.5, 1.0, 1.0]
    prior_lower = [0.5, 0.0, -1.0]
    stepsize = [0.007, 0.05, 0.05]

    # Initialize parameters
    parameters = Float64[Random.rand(Distributions.Uniform(i,k)) for (i,k) in zip(prior_lower,prior_upper)]
    new_parameters = []

    chi2 = likelihood(sn_data, parameters)
    new_chi2 = Inf

    best_chi2 = Inf
    best_param = []
    parameters_trace = hcat([chi2, parameters[1], parameters[2], parameters[3]])

    for i in range(1, max_step)

        new_parameters = parameters .+ Random.randn(3) .* stepsize

        if all(prior_lower.<new_parameters.<prior_upper)
            new_chi2 = likelihood(sn_data, new_parameters)
        else
            new_chi2 = Inf
        end

        if new_chi2 < chi2
            chi2 = new_chi2
            parameters = new_parameters
            parameters_trace = hcat(parameters_trace, [chi2, parameters[1], parameters[2], parameters[3]])
        elseif exp(-(new_chi2-chi2)/2.0) > Random.rand()
            chi2 = new_chi2
            parameters = new_parameters
            parameters_trace = hcat(parameters_trace, [chi2, parameters[1], parameters[2], parameters[3]])
        else
            nothing
        end

        if new_chi2 < best_chi2
            best_chi2 = chi2
            best_param = parameters
        end

    end

    return best_param, best_chi2, parameters_trace
end