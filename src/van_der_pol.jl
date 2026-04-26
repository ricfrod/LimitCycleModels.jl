"""
van der Pol, Balth. (1926). LXXXVIII. On “relaxation-oscillations.” The London, Edinburgh,
and Dublin Philosophical Magazine and Journal of Science, 2(11), 978–992.
https://doi.org/10.1080/14786442608564127
"""

@kwdef struct VanDerPol <: AbstractLimitCycleModel
    μ::Float64
end

function (self::VanDerPol)(
        du::AbstractVector{Float64},
        u::AbstractVector{Float64},
        args...
    )::Nothing

    x::Float64, y::Float64 = u
    (; μ) = self

    du[1] = y
    du[2] = μ * (1 - x^2) * y - x

    return nothing
end
