let
    du = zeros(8)
    fowm = fowm! = FowmFriction()
    fowm!(du, fowm.u0, [fowm.ω_gc, fowm.z], nothing)

    expected_du = [
        0.2431409645038534
        0.18842703701689378
        5.653732284412239
        0.06242929005890585
        0.8538298909656432
        7.895963436419567
        -1.0078760101441913
        -392.14654970585804
    ]

    @test all(isapprox.(du, expected_du))
end
