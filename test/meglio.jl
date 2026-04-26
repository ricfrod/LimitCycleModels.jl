let
    du = zeros(3)
    p = Meglio()
    u0 = [6088.499179683961, 282.4512708285756, 5089.462041088525]

    model!(du, u0, p, nothing)

    expected_du = [
        0.18039999999999998
        0.006254745728489275
        0.33781332063915137
    ]

    @test all(isapprox.(du, expected_du))
end
