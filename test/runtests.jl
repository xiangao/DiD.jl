using Test
using ETWFE
using DataFrames

@testset "DiD.jl (deprecated)" begin
    # This package is deprecated in favor of Panelest.jl's etwfe()+emfx().
    # att_gt never estimated anything — it returned hardcoded placeholder
    # numbers. It now errors instead of silently returning fake data; this
    # is the only behavior left worth testing.
    df = DataFrame(
        id    = repeat(1:30, inner = 3),
        time  = repeat(1:3, outer = 30),
        group = vcat(repeat([0], 10 * 3), repeat([2], 10 * 3), repeat([3], 10 * 3)),
        y     = randn(90),
    )

    @test_throws ErrorException att_gt(df, :y, :time, :id, :group)
end
