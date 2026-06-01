module ManifoldsBaseMakiePlotsExt
#
#
# For the case where one actually loads both plotting libraries,
# this uses Preferences.jl to persistently specify which backend to use
using ManifoldsBase
import ManifoldsBase: splot_slope, plot_check_geodesic

function _preferred_plotting_ext()
    if Symbol(ManifoldsBase.get_plotting_backend()) === :Makie
        return Base.get_extension(ManifoldsBase, :ManifoldsBaseMakieExt)
    else
        return Base.get_extension(ManifoldsBase, :ManifoldsBasePlotsExt)
    end
end
function ManifoldsBase.plot_slope(
        x, y; slope = 2, line_base = 0, a = 0, b = 2.0, i = 1, j = length(x),
    )
    return _preferred_plotting_ext()._plot_slope(x, y; slope, line_base, a, b, i, j)
end
function ManifoldsBase.plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
    return _preferred_plotting_ext()._plot_check_geodesic(T, N, e_norm, e_pt, e_alpha)
end

end # module ManifoldsBaseMakiePlotsExt
