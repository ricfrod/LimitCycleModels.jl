"""
Default parameters from:
Rodrigues, R.F., Trierweiler, J.O. & Farenzena, M. Parameter estimation and statistical
analysis of limit cycle models using a local optimization algorithm: applications for
simplified well models of slugging oil systems. Braz. J. Chem. Eng. (2026).
https://doi.org/10.1007/s43153-025-00628-2
"""
@kwdef struct Fowm <: AbstractLimitCycleModel
    ωᵤ::Real = 2.802
    m_lstill::Real = 6.442e1
    C_g::Real = 1.12e-3
    ϵ::Real = 6.617e-1
    V_eb::Real = 6.105e1
    C_out::Real = 2.039e-3
    K_w::Real = 6.876e-4
    Kₐ::Real = 2.293e-5
    Kᵣ::Real = 1.269e2
    ρ_l::Real = 900.0
    T::Real = 298.0
    M::Real = 18.0
    θ::Real = π / 4
    Pₛ::Real = 1013250.0
    Pᵣ::Real = 2.25e7
    α_gw::Real = 0.0188
    ρ_mres::Real = 891.9523
    Lᵣ::Real = 4497.0
    Lₜ::Real = 1639.0
    Lₐ::Real = 1118.0
    Dᵣ::Real = 0.152
    Dₜ::Real = 0.15
    Dₐ::Real = 0.14
    Hₜ::Real = 1279.0
    H_pdg::Real = 1117.0
    H_vgl::Real = 916.0
    z::Real = 0.84
    ω_gc::Real = 1.27
    R::Real = 8314.46261815324
    g::Real = 9.80665
end


function get_model(p::Fowm)
    (
        ; ρ_l, T, M, θ, Pₛ, Pᵣ, α_gw,
        ρ_mres, Lᵣ, Lₜ, Lₐ, Dᵣ, Dₜ,
        Dₐ, Hₜ, H_pdg, H_vgl, ωᵤ,
        m_lstill, C_g, C_out, V_eb,
        ϵ, K_w, Kₐ, Kᵣ, R, g,
    ) = p

    A = Dᵣ^2 * π / 4
    Vᵣ = ωᵤ * Lᵣ * A
    Vₜ = Lₜ * Dₜ^2 * π / 4
    Vₐ = Lₐ * Dₐ^2 * π / 4
    α_lw = 1 - α_gw

    function model!(du, u, p, _)

        m_gb, m_gr, m_lr, m_ga, m_gt, m_lt = u
        ω_gc, z = p

        # Pressure
        ρ_gr = m_gr / (Vᵣ - (m_lr + m_lstill) / ρ_l)
        V_gt = Vₜ - m_lt / ρ_l
        ρ_gt = m_gt / V_gt
        ρ̄ₜ = (m_gt + m_lt) / Vₜ

        P_eb = m_gb * R * T / (M * V_eb)
        P_rt = ρ_gr * R * T / M
        P_rb = P_rt + (m_lr + m_lstill) * g * sin(θ) / A
        P_tt = (ρ_gt * R * T) / M
        P_tb = P_tt + ρ̄ₜ * g * H_vgl
        P_pdg = P_tb + ρ_mres * g * (H_pdg - H_vgl)
        P_bh = P_pdg + ρ_mres * g * (Hₜ - H_pdg)
        P_ai = (R * T / (Vₐ * M) + (g * Lₐ / Vₐ)) * m_ga

        # Mass flow
        α_gr = m_gr / (m_gr + m_lr)
        α_lr = 1 - α_gr
        α_gt = m_gt / (m_gt + m_lt)
        α_lt = 1 - α_gt
        ρ_ai = M * P_ai / (R * T)
        ω_out = C_out * z * sqrt(max(0.0, ρ_l * ((P_rt - Pₛ))))
        ω_wh = K_w * sqrt(max(0.0, ρ_l * ((P_tt - P_rb))))
        ωᵣ = Kᵣ * max(0.0, 1 - 0.2 * P_bh / Pᵣ - 0.8 * (P_bh / Pᵣ)^2)

        ω_gout = α_gr * ω_out
        ω_lout = α_lr * ω_out
        ω_gwh = α_gt * ω_wh
        ω_lwh = α_lt * ω_wh
        ω_g = C_g * max(0.0, P_eb - P_rb)
        ω_iv = Kₐ * sqrt(max(0.0, ρ_ai * (P_ai - P_tb)))
        ω_gin = α_gw * ωᵣ
        ω_lin = α_lw * ωᵣ

        # Derivative
        du[1] = (1 - ϵ) * (ω_gwh) - ω_g
        du[2] = ϵ * (ω_gwh) + ω_g - ω_gout
        du[3] = ω_lwh - ω_lout
        du[4] = ω_gc - ω_iv
        du[5] = ω_gin + ω_iv - ω_gwh
        du[6] = ω_lin - ω_lwh

        return nothing
    end

    return model!
end
