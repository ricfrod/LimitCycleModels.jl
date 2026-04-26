"""
Lotka, A. J. (1920). Analytical Note on Certain Rhythmic Relations in Organic Systems.
Proceedings of the National Academy of Sciences, 6(7), 410–415.
https://doi.org/10.1073/pnas.6.7.410

Volterra, V. (1928). Variations and Fluctuations of the Number of Individuals in Animal
Species living together. ICES Journal of Marine Science, 3(1), 3–51.
https://doi.org/10.1093/icesjms/3.1.3
"""

@kwdef struct LotkaVolterra <: AbstractLimitCycleModel
    α::Float64
    β::Float64
    δ::Float64
    γ::Float64
end


function (self::LotkaVolterra)(
        du::AbstractVector{Float64},
        u::AbstractVector{Float64},
        args...
    )::Nothing

    x, y = u
    (; α, β, δ, γ) = self

    du[1] = α * x - β * x * y
    du[2] = δ * x * y - γ * y

    return nothing
end
