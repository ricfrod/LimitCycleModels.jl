@kwdef struct FowmFriction <: AbstractLimitCycleModel
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
    A::Real = 0.0
    Vᵣ::Real = 0.0
    Vₜ::Real = 0.0
    Vₐ::Real = 0.0
    α_lw::Real = 0.0
    ε::Real = 1.0e-7
    μ_l::Real = 1.43e-4
    μ_g::Real = 1.39e-5
    ξ::Float64 = eps(Float64)
    u0::AbstractVector{<:Real} = [
        4890.431101521542
        1470.428382621798
        28158.468433445065
        3479.00609623567
        1013.4944459459704
        15726.74205428866
        0.0
        0.0
    ]

    function FowmFriction(
            ωᵤ,
            m_lstill,
            C_g,
            ϵ,
            V_eb,
            C_out,
            K_w,
            Kₐ,
            Kᵣ,
            ρ_l,
            T,
            M,
            θ,
            Pₛ,
            Pᵣ,
            α_gw,
            ρ_mres,
            Lᵣ,
            Lₜ,
            Lₐ,
            Dᵣ,
            Dₜ,
            Dₐ,
            Hₜ,
            H_pdg,
            H_vgl,
            z,
            ω_gc,
            R,
            g,
            A,
            Vᵣ,
            Vₜ,
            Vₐ,
            α_lw,
            ε,
            μ_l,
            μ_g,
            ξ,
            u0,
        )
        A = Dᵣ^2 * π / 4
        Vᵣ = ωᵤ * Lᵣ * A
        Vₜ = Lₜ * Dₜ^2 * π / 4
        Vₐ = Lₐ * Dₐ^2 * π / 4
        α_lw = 1 - α_gw
        return new(
            ωᵤ,
            m_lstill,
            C_g,
            ϵ,
            V_eb,
            C_out,
            K_w,
            Kₐ,
            Kᵣ,
            ρ_l,
            T,
            M,
            θ,
            Pₛ,
            Pᵣ,
            α_gw,
            ρ_mres,
            Lᵣ,
            Lₜ,
            Lₐ,
            Dᵣ,
            Dₜ,
            Dₐ,
            Hₜ,
            H_pdg,
            H_vgl,
            z,
            ω_gc,
            R,
            g,
            A,
            Vᵣ,
            Vₜ,
            Vₐ,
            α_lw,
            ε,
            μ_l,
            μ_g,
            ξ,
            u0,
        )
    end
end


function (self::FowmFriction)(
        du::AbstractVector{<:Real},
        u::AbstractVector{<:Real},
        p::AbstractVector{<:Real},
        _
    )::Nothing

    (
        ; ρ_l, T, M, θ, Pₛ, Pᵣ, α_gw,
        ρ_mres, Lᵣ, Lₜ, Lₐ, Dᵣ, Dₜ,
        Hₜ, H_pdg, H_vgl, α_lw,
        m_lstill, C_g, C_out, V_eb,
        ϵ, K_w, Kₐ, Kᵣ, ε, μ_l, μ_g,
        A, Vᵣ, Vₜ, Vₐ, R, g, ξ
    ) = self
    m_gb, m_gr, m_lr, m_ga, m_gt, m_lt, ΔPₜ, ΔPᵣ = u
    ω_gc, z = p

    # Pressure
    ρ_gr = m_gr / (Vᵣ - (m_lr + m_lstill) / ρ_l)
    V_gt = Vₜ - m_lt / ρ_l
    ρ_gt = m_gt / V_gt
    ρ̄ₜ = (m_gt + m_lt) / Vₜ

    P_eb = m_gb * R * T / (M * V_eb)
    P_rt = ρ_gr * R * T / M
    P_rb = P_rt + ΔPᵣ + (m_lr + m_lstill) * g * sin(θ) / A
    P_tt = (ρ_gt * R * T) / M
    P_tb = P_tt + ΔPₜ + ρ̄ₜ * g * H_vgl
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

    # Friction loss
    u_sgt = 4 * (ω_iv + ω_gin) / (ρ_gt * π * Dₜ^2)
    u_slt = 4 * ω_lin / (ρ_l * π * Dₜ^2)
    ūₜ = u_sgt + u_slt
    μ̄ₜ = α_gt * μ_g + α_lt * μ_l
    Reₜ = abs(ρ̄ₜ * ūₜ * Dₜ / μ̄ₜ)
    inner_logt = max(ξ, ((ε / Dₜ) / 3.7)^1.11 + 6.9 / (Reₜ + ξ))
    fₜ = (-1.8 * log10(inner_logt))^(-2)

    u_sgr = 4 * (ω_iv + ω_gin) / (ρ_gr * π * Dᵣ^2)
    u_slr = 4 * ω_lin / (ρ_l * π * Dᵣ^2)
    ūᵣ = u_sgr + u_slr
    μ̄ᵣ = α_gr * μ_g + α_lr * μ_l
    ρ̄ᵣ = (m_gr + m_lr) / Vᵣ
    Reᵣ = abs(ρ̄ᵣ * ūᵣ * Dᵣ / μ̄ᵣ)
    inner_logr = max(ξ, ((ε / Dᵣ) / 3.7)^1.11 + 6.9 / (Reᵣ + ξ))
    fᵣ = (-1.8 * log10(inner_logr))^(-2)

    # Derivatives
    du[1] = (1 - ϵ) * (ω_gwh) - ω_g
    du[2] = ϵ * (ω_gwh) + ω_g - ω_gout
    du[3] = ω_lwh - ω_lout
    du[4] = ω_gc - ω_iv
    du[5] = ω_gin + ω_iv - ω_gwh
    du[6] = ω_lin - ω_lwh
    du[7] = ΔPₜ - fₜ * Lₜ * ūₜ^2 / Dₜ / ρ̄ₜ
    du[8] = ΔPᵣ - fᵣ * Lᵣ * ūᵣ^2 / Dᵣ / ρ̄ᵣ

    return nothing
end
