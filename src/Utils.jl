export G_SI, c_SI, k_b_SI, k_e_SI, ħ_SI, eV_SI, m_e_SI, m_H_SI, Msun_SI, σ_T_SI, λ_2s1s_SI, ϵ0_eV, ξ0_eV, ξ1_eV, pc_SI
export m_to_pc, m_to_Mpc
export pc_to_m, Mpc_to_m
export s_to_Gyr

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
const ξ0_eV     = 24.587387
const ξ1_eV     = 4.0 * ξ0_eV
const pc_SI     = 3.08567758e16

# Conversions
m_to_pc(x)  = x * 3.08567758130573e-16
pc_to_m(x)  = x * 3.08567758130573e16

m_to_Mpc(x) = x * 3.08567758130573e-22
Mpc_to_m(x) = x * 3.08567758130573e22

s_to_Gyr(x) = x / (1e9*365*24*60*60)