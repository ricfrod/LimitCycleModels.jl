"""
Default parameters from:
di Meglio, F., Alstad, V., Petit, N., & Kaasa, G.-O. (2012). Stabilization of slugging in
oil production facilities ωith or ωithout upstream pressure sensors. Journal of Process
Control, 22(4), 809–822. https://doi.org/10.1016/j.jprocont.2012.02.014
"""
@kwdef struct Meglio <: AbstractLimitCycleModel
    z::Real = 0.3
    m_lstill::Real = 3.73e4
    C_g::Real = 1.0e-4
    ϵ::Real = 0.78
    V_eb::Real = 48.0
    C_out::Real = 2.8e-3
    ρₗ::Real = 900.0
    T::Real = 363.0
    M::Real = 22.0
    θᵣ::Real = π / 4
    Pₛ::Real = 6.6e5
    A::Real = pi * (1.5e-1 / 2)^2 #4.4e-3
    Vᵣ::Real = 5200 * pi * (1.5e-1 / 2)^2
    ω_ing::Real = 8.2e-1
    ω_inl::Real = 11.75
    R::Real = 8314.46261815324
    g::Real = 9.80665
end


function model!(
        du::AbstractVector{Float64},
        u::AbstractVector{Float64},
        p::Meglio,
        _,
    )::Nothing
    m_eb, m_rg, m_rl = u

    (;
        z, m_lstill, C_g, ϵ, V_eb, C_out, ρₗ, T, M, θᵣ,
        Pₛ, A, Vᵣ, ω_ing, ω_inl, R, g,
    ) = p

    P_eb = m_eb * R * T / (M * V_eb)
    P_rt = m_rg * R * T / (M * (Vᵣ - (m_rl + m_lstill) / ρₗ))
    P_rb = P_rt + (m_rl + m_lstill) * g * sin(θᵣ) / A

    ω_g = C_g * max(0, P_eb - P_rb)
    ω_out = z * C_out * sqrt(ρₗ * max(0, P_rt - Pₛ))
    ω_outg = (m_rg / m_rl) * ω_out

    du[1] = (1 - ϵ) * ω_ing - ω_g
    du[2] = ϵ * ω_ing + ω_g - ω_outg
    du[3] = ω_inl - ω_out

    return nothing
end
