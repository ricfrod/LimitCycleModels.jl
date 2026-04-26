"""
Default parameters from:
Apio, A., Botelho, V. R., & Trierweiler, J. O. (2018). Parameter estimation of models with
limit cycle based on the reformulation of the objective function. Computers & Chemical
Engineering, 109, 236–248. https://doi.org/10.1016/j.compchemeng.2017.11.009
"""
@kwdef struct Jobses
    C_eo::Real = 0
    C_Po::Real = 0
    C_xo::Real = 0
    C_So::Real = 200
    D_f::Real = 0.07
    K_S::Real = 0.5
    k_1::Real = 16
    k_2::Real = 0.497
    k_3::Real = 0.00383
    μ_max::Real = 1
    m_S::Real = 2.16
    m_P::Real = 1.3
    Y_Px::Real = 0.06
    Y_Sx::Real = 0.02339905
end

function (self::Jobses)(
    du::AbstractVector{Float64},
    u::AbstractVector{Float64},
    args...
)::Nothing

    (
        ; C_eo, C_Po, C_xo, C_So, D_f, K_S, k_1, k_2,
        k_3, μ_max, m_S, m_P, Y_Px, Y_Sx,
    ) = self

    C_e, C_P, C_x, C_S = u

    _x = (C_S * C_e / (K_S + C_S))

    du[1] = (
        μ_max * _x * (-1 / Y_Sx)
            - m_S * C_x
            + D_f * (C_So - C_S)
    )
    du[2] = (
        μ_max * _x
            + D_f * (C_xo - C_x)
    )
    du[3] = (
        μ_max * _x * (k_1 - k_2 * C_P + k_3 * C_P^2)
            + D_f * (C_eo - C_e)
    )
    du[4] = (
        μ_max * _x * (1 / Y_Px)
            + m_P * C_x
            + D_f * (C_Po - C_P)
    )

    return nothing
end
