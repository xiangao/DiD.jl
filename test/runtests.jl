using Test
using ETWFE
using DataFrames

@testset "DiD.jl" begin
    # Create simple panel data
    # 3 periods: 1, 2, 3
    # Group 0: Never treated
    # Group 2: Treated in period 2
    # Group 3: Treated in period 3
    
    df = DataFrame(
        id = repeat(1:30, inner=3),
        time = repeat(1:3, outer=30),
        group = vcat(repeat([0], 10*3), repeat([2], 10*3), repeat([3], 10*3)),
    )
    
    # Baseline
    df.y = randn(90)
    
    # Treatment effects
    # Group 2, time >= 2: +5.0
    # Group 3, time >= 3: +10.0
    
    for i in 1:nrow(df)
        if df.group[i] == 2 && df.time[i] >= 2
            df.y[i] += 5.0
        elseif df.group[i] == 3 && df.time[i] >= 3
            df.y[i] += 10.0
        end
    end
    
    res = att_gt(df, :y, :time, :id, :group)
    
    @test length(res.group) > 0
    @test length(res.att) > 0
    
    # event-study
    es = emfx(res, type="event")
    @test "event" in names(es)
    
    # calendar
    cal = emfx(res, type="calendar")
    @test "calendar_time" in names(cal)
    
    # overall
    ov = emfx(res, type="overall")
    @test "estimate" in names(ov)
end
