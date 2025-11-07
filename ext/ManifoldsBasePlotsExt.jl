module ManifoldsBasePlotsExt

using ManifoldsBase
using Plots
using Printf: @sprintf
import ManifoldsBase: plot_slope

function ManifoldsBase.plot_slope(
        x,
        y;
        slope = 2,
        line_base = 0,
        a = 0,
        b = 2.0,
        i = 1,
        j = length(x),
    )
    fig = plot(
        x,
        y;
        xaxis = :log,
        yaxis = :log,
        label = "\$E(t)\$",
        linewidth = 3,
        legend = :topleft,
        color = :lightblue,
    )
    s_line = [exp10(line_base + t * slope) for t in log10.(x)]
    plot!(
        fig,
        x,
        s_line;
        label = "slope s=$slope",
        linestyle = :dash,
        color = :black,
        linewidth = 2,
    )
    if (i != 0) && (j != 0)
        best_line = [exp10(a + t * b) for t in log10.(x[i:j])]
        plot!(
            fig,
            x[i:j],
            best_line;
            label = "best slope $(@sprintf("%.4f", b))",
            color = :blue,
            linestyle = :dot,
            linewidth = 2,
        )
    end
    return fig
end
function ManifoldsBase.plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
    fig = plot(
        T[1:N-1], e_norm;
        label = "\$\\|X_i\\|-c\$", linewidth = 3, legend = :topright, color = :lightblue,
    )
    plot!(
        fig, T[1:N-2], e_pt; label = "\$E_{\\mathrm{n}})\$", linewidth = 3, color = :blue,
    )
    plot!(
        fig, T[1:N-2], e_alpha; label = "\$E_α\$", linewidth = 3, color = :green,
    )
    return fig
end
end
