let
    du = zeros(3)
    f! = Meglio()
    u0 = [6088.499179683961, 282.4512708285756, 5089.462041088525]

    f!(du, u0, nothing, nothing)

    expected_du = [
        0.18039999999999998
        0.006254745728489275
        0.33781332063915137
    ]

    @test all(isapprox.(du, expected_du))
end
