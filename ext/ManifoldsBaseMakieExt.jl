module ManifoldsBaseMakieExt

using ManifoldsBase
using Makie
using Printf: @sprintf
import ManifoldsBase: splot_slope, plot_check_geodesic

function ManifoldsBase.plot_slope(
        x, y; slope = 2, line_base = 0, a = 0, b = 2.0, i = 1, j = length(x),
    )
    return _plot_slope(x, y; slope = slope, line_base = line_base, a = a, b = b, i = i, j = j)
end
function _plot_slope(x, y; slope = 2, line_base = 0, a = 0, b = 2.0, i = 1, j = length(x))
    fig = Makie.Figure()
    # Setup the log log plot
    ax = Makie.Axis(
        fig[1, 1];
        xscale = log10, yscale = log10, xlabel = "t", ylabel = L"E(t)", title = "Slope plot",
    )
    # Main error data
    Makie.lines!(ax, x, y; label = L"E(t)", linewidth = 3, color = :lightblue)
    # Reference slope
    s_line = [exp10(line_base + t * slope) for t in log10.(x)]
    Makie.lines!(
        ax, x, s_line; label = "slope s=$slope",
        linestyle = :dash, color = :black, linewidth = 2
    )

    # Best-fit slope segment if i and j are provided
    if (i != 0) && (j != 0)
        best_line = [exp10(a + t * b) for t in log10.(x[i:j])]
        Makie.lines!(
            ax, x[i:j], best_line;
            label = "best slope $(@sprintf("%.4f", b))",
            color = :blue, linestyle = :dot, linewidth = 2,
        )
    end
    Makie.axislegend(ax; position = :lt)
    return fig
end
function ManifoldsBase.plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
    return _plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
end
function _plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xlabel = "t", ylabel = "error magnitude", title = "Geodesic checks",
    )
    # Speed deviation: |‖Xᵢ‖ − mean‖X‖|
    Makie.lines!(
        ax, T[1:(N - 1)], e_norm;
        label = L"|\|X_i\| - \overline{\|X\|}|", linewidth = 3, color = :lightblue
    )
    # Parallel-transport mismatch (norm): ‖Xᵢ − PT_{i+1→i}(X_{i+1})‖
    lines!(
        ax, T[1:(N - 2)], e_pt;
        label = "PT difference (norm)", linewidth = 3, color = :blue
    )
    # Parallel-transport mismatch (angle): |1 − cos∠(Xᵢ, PT_{i+1→i}X_{i+1})|
    lines!(
        ax, T[1:(N - 2)], e_alpha;
        label = "PT difference (angle)", linewidth = 3, color = :green
    )
    axislegend(ax; position = :rt)
    return fig
end
end
