using LimitCycleModels
using Test

@testset "LimitCycleModels.jl" begin

    @testset "van der Pol" begin
        include("van_der_pol.jl")
    end

    @testset "Lotka-Volterra" begin
        include("lotka_volterra.jl")
    end

    @testset "Jobses" begin
        include("jobses.jl")
    end

    @testset "Meglio" begin
        include("meglio.jl")
    end

    @testset "Fowm" begin
        include("fowm.jl")
    end

    @testset "Fowm Friction" begin
        include("fowm_friction.jl")
    end

end
