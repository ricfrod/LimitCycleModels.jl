let
    p = LotkaVolterra(; α = 1.5, β = 1.0, δ = 1.0, γ = 3.0)
    du = zeros(2)
    u = [1.0, 1.0]
    model!(du, u, p, 0.0)

    @test all(isapprox.(du, [0.5, -2.0]))
end
