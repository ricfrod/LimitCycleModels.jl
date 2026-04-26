let
    f! = Jobses(;
        C_eo = 0,
        C_Po = 0,
        C_xo = 0,
        C_So = 200,
        D_f = 0.07,
        K_S = 0.5,
        k_1 = 16,
        k_2 = 0.497,
        k_3 = 0.00383,
        μ_max = 1,
        m_S = 2.16,
        m_P = 1.3,
        Y_Px = 0.06,
        Y_Sx = 0.02339905,
    )

    du = zeros(4)
    u = zeros(4)
    f!(du, u, nothing, nothing)

    @test all(isapprox.(du, [14, 0, 0, 0]))

end
