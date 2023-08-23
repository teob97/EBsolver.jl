export RecombinationHistory
export  Xe_saha_equation

Base.@kwdef struct RecombinationHistory

    cosmo :: BackgroundCosmology = BackgroundCosmology()

    # Helium fraction
    Yp :: Float64 = 0.245

    # Start and end points for recombination arrays
    x_start :: Float64 = cosmo.x_start
    x_end :: Float64 = cosmo.x_end
    
    # Numbers of points of Xe,ne array
    npts_rec_arrays :: Int64 = 4000

    # Limit for when to switch between Saha and Peebles
    Xe_saha_limit :: Float64 = 0.99

end

function Xe_saha_equation(RH::RecombinationHistory, x::Float64)
    n_b :: Float64 = RH.cosmo.Ω0_B * eval_ρ0_crit(RH.cosmo.H0_SI) * exp(-3x) / m_H_SI
    T_b :: Float64 = RH.cosmo.T0_CMB * exp(-x)
    k :: Float64 = (ħ_SI^3 * n_b)^(-1) * (m_e_SI*k_b_SI*T_b/2π)^(3/2) * exp(-ϵ0_SI/(k_b_SI*T_b))
    # Avoid subtract huge floats with same magnitude 
    if 4/k < 1e-5
        Xe = 1.0
    else 
        Xe = 0.5 * (sqrt(k^2+4k) - k)
    end 
    return Xe
end