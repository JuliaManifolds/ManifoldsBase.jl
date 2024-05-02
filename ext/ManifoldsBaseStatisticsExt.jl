module ManifoldsBaseStatisticsExt

if isdefined(Base, :get_extension)
    using Statistics
    using ManifoldsBase
    import ManifoldsBase: find_best_slope_window
else
    # imports need to be relative for Requires.jl-based workflows:
    # https://github.com/JuliaArrays/ArrayInterface.jl/pull/387
    using ..Statistics
    using ..ManifoldsBase
    import ..ManifoldsBase: find_best_slope_window
end

function ManifoldsBase.find_best_slope_window(
    X,
    Y,
    window = nothing;
    slope::Real = 2.0,
    slope_tol::Real = 0.1,
)
    n = length(X)
    if window !== nothing && (any(window .> n))
        error(
            "One of the window sizes ($(window)) is larger than the length of the signal (n=$n).",
        )
    end
    a_best = 0
    b_best = -Inf
    i_best = 0
    j_best = 0
    r_best = 0 # longest interval
    for w in (window === nothing ? (2:n) : [window...])
        for j in 1:(n - w + 1)
            x = X[j:(j + w - 1)]
            y = Y[j:(j + w - 1)]
            # fit a line a + bx
            c = cor(x, y)
            b = std(y) / std(x) * c
            a = mean(y) - b * mean(x)
            # look for the largest interval where b is within slope tolerance
            r = (maximum(x) - minimum(x))
            if (r > r_best) && abs(b - slope) < slope_tol #longer interval found.
                r_best = r
                a_best = a
                b_best = b
                i_best = j
                j_best = j + w - 1 #last index (see x and y from before)
            end
            # not best interval - maybe it is still the (first) best slope?
            if r_best == 0 && abs(b - slope) < abs(b_best - slope)
                # but do not update `r` since this indicates only a best r
                a_best = a
                b_best = b
                i_best = j
                j_best = j + w - 1 #last index (see x and y from before)
            end
        end
    end
    return (a_best, b_best, i_best, j_best)
end
end
