let
    p = VanDerPol(1.0)
    du = zeros(2)
    u = [2.0, 0.0]
    model!(du, u, p, 0.0)

    @test all(isapprox.(du, [0.0, -2.0]))
end
