module ManifoldsBasePlotsExt

using ManifoldsBase
using Plots
using Printf: @sprintf

function ManifoldsBase.plot_slope(
        ::Val{:Plots}, x, y;
        slope = 2, line_base = 0, a = 0, b = 2.0, i = 1, j = length(x),
    )
    fig = plot(
        x, y; xaxis = :log, yaxis = :log, label = "\$E(t)\$",
        linewidth = 3, legend = :topleft, color = :lightblue,
    )
    s_line = [exp10(line_base + t * slope) for t in log10.(x)]
    plot!(
        fig, x, s_line; label = "slope s=$slope",
        linestyle = :dash, color = :black, linewidth = 2,
    )
    if (i != 0) && (j != 0)
        best_line = [exp10(a + t * b) for t in log10.(x[i:j])]
        plot!(
            fig, x[i:j], best_line;
            label = "best slope $(@sprintf("%.4f", b))",
            color = :blue, linestyle = :dot, linewidth = 2,
        )
    end
    return fig
end
function ManifoldsBase.plot_check_geodesic(::Val{:Plots}, T, N, e_norm, e_pt, e_alpha)
    fig = plot(
        T[1:(N - 1)], e_norm;
        # Deviation from constant speed: |‖Xᵢ‖ − mean‖X‖|
        label = "Speed deviation: |‖Xᵢ‖ − mean‖X‖|",
        linewidth = 3, legend = :topright, color = :lightblue,
        xlabel = "t", ylabel = "error magnitude", title = "Geodesic checks"
    )
    plot!(
        fig, T[1:(N - 2)], e_pt;
        # Parallel transport mismatch: ‖Xᵢ − PT_{i+1→i}(X_{i+1})‖
        label = "PT difference (norm)",
        linewidth = 3, color = :blue,
    )
    plot!(
        fig, T[1:(N - 2)], e_alpha;
        # Angle-preservation error: |1 − cos∠(Xᵢ, PT_{i+1→i}X_{i+1})|
        label = "PT difference (angle)",
        linewidth = 3, color = :green,
    )
    return fig
end
# Check whether the default fallback is not Plots and if so set it to this.
# Then, if no persistent backend is set, plotting “just works”
ManifoldsBase._MANIFOLDSBASE_PLOTTING_BACKEND != "Plots" && (ManifoldsBase.set_plotting_backend!("Plots"; only_fallback = true))
nm = isnothing(Base.get_extension(ManifoldsBase, :ManifoldsBaseMakieExt))
pb = ManifoldsBase.get_plotting_backend()
# Makie is not loaded but current default – warn
(pb != "Plots") && nm && (@warn "PLots was loaded, but the current backend is `Makie`, which is not loaded. consider calling `ManifoldsBase.set_plotting_backend!(\"Plots\")`.")
end
