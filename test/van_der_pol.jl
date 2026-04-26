let
    f! = VanDerPol(1.0)
    du = zeros(2)
    u = [2.0, 0.0]
    f!(du, u, nothing, nothing)

    @test all(isapprox.(du, [0.0, -2.0]))
end
