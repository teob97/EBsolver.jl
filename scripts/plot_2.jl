using EBsolver
using Plots
using Logging

let

    @info "Setting up Background Cosmology"

    path = "scripts/plots/RecombinationHistory/milestone2.png"
    rh = RecombinationHistory()
    bc = BackgroundCosmology()
    x = range(rh.cosmo.x_start, -1e-2, 4000)

    @info "Making τ(x) and derivatives graph"

    p1 = plot(
        x, 
        τ_of_x(rh, false).(x),
        yaxis=:log,
        xlim=[rh.x_start, rh.x_end],
        ylim=[10^-8, 10^8],
        linewidth=2,
        xlabel="x=log(a)",
        label = "τ(x)"
    )

    p1 = plot!(
        x, 
        -dτdx_of_x(rh, false).(x),
        linewidth=2,
        xlabel="x=log(a)",
        label = "-τ'(x)"
    )

    p1 = plot!(
        x, 
        ddτddx_of_x(rh, false).(x),
        linewidth=2,
        xlabel="x=log(a)",
        label = "τ''(x)"
    )

    @info "Making τ(x) and derivatives graph with reionization"

    p2 = plot(
        x, 
        τ_of_x(rh, true).(x),
        yaxis=:log,
        xlim=[rh.x_start, rh.x_end],
        ylim=[10^-8, 10^8],
        linewidth=2,
        xlabel="x=log(a)",
        label = "τ(x)"
    )

    p2 = plot!(
        x, 
        -dτdx_of_x(rh, true).(x),
        linewidth=2,
        xlabel="x=log(a)",
        label = "-τ'(x)"
    )

    p2 = plot!(
        x, 
        ddτddx_of_x(rh, true).(x),
        linewidth=2,
        xlabel="x=log(a)",
        label = "τ''(x)"
    )

    @info "Making X_e(x) graph"

    p3 = plot(
        x,
        Xe_of_x(rh, true).(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="X_e(x) with reionization",
        line=:dash,
    )

    p3 = plot!(
        x,
        Xe_of_x(rh, false).(x),
        linewidth=2,
        label="X_e(x)"
        )

    @info "Making visibility function g(x) graph"
    
    g_noreion = visibility_function_of_x(rh, false)
    dg_noreion = diff(g_noreion)
    ddg_noreion = diff(dg_noreion)

    g = visibility_function_of_x(rh, true)
    dg = diff(g)
    ddg = diff(dg)

    p4 = plot(
        x,
        g_noreion.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g(x)",
        line=:dash,
    )

    p4 = plot!(
        x,
        g.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g(x) with reionization",
    )
    
    @info "Making g'(x) graph"

    p5 = plot(
        x,
        dg_noreion.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g'(x)",
        line=:dash,
    )

    p5 = plot!(
        x,
        dg.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g'(x) with reionization",
    )
    
    @info "Making g''(x) graph"

    p6 = plot(
        x,
        ddg_noreion.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g''(x)",
        line=:dash,
    )

    p6 = plot!(
        x,
        ddg.(x),
        xlim=[rh.x_start, rh.x_end],
        linewidth=2,
        xlabel="x=log(a)",
        label="g''(x) with reionization",
    )

    plot(p1, p2, p3, p4, p5, p6, size = (1000, 600))
    savefig(path)

    @info "Image saved at "*path

end