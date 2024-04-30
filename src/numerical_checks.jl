
@doc raw"""
    check_retraction(
        M::AbstractManifold,
        rectraction_method::AbstractRetractionMethod,
        p=rand(M),
        X=rand(M; vector_at=p);
        #
        exactness_tol = 1e-12,
        io = nothing,
        limits = (-8.0,1),
        log_range = range(limits[1], limits[2]; length=N),
        N = 101,
        name = "retraction",
        plot = false,
        second_order = true
        slope_tol = 0.1,
        error = :none,
        window = nothing,
    )

Check numerically wether the retraction `retraction_method` is correct.
This requires the [`exp`](@ref) function to be implemented for the [`AbstractManifold`](@ref) `M`.

This implements a method similar to [Boumal:2023; Section 4.8 or Section 6.8](@cite).

Note that if the errors are below the given tolerance and the method is exact,
no plot is generated,

# Keyword arguments

* `exactness_tol`:     if all errors are below this tolerance, the differential is considered to be exact
* `io`:                provide an `IO` to print the result to
* `limits`:            specify the limits in the `log_range`, that is the exponent for the range
* `log_range`:         specify the range of points (in log scale) to sample the differential line
* `N`:                 number of points to verify within the `log_range` default range ``[10^{-8},10^{0}]``
* `name`:              name to display in the plot
* `plot`:              whether to plot the result (if `Plots.jl` is loaded).
  The plot is in log-log-scale. This is returned and can then also be saved.
* `second_order`:      check whether the retraction is of second order. if set to `false`, first order is checked.
* `slope_tol`:         tolerance for the slope (global) of the approximation
* `error`:             specify how to report errors: `:none`, `:info`, `:warn`, or `:error` are availablem, cf. [is_point](@ref)
* `window`:            specify window sizes within the `log_range` that are used for the slope estimation.
  the default is, to use all window sizes `2:N`.
"""
function check_retraction(
    M::AbstractManifold,
    retraction_method::AbstractRetractionMethod,
    p = rand(M),
    X = rand(M; vector_at = p);
    exactness_tol = 1e-12,
    io::Union{IO,Nothing} = nothing,
    limits = (-8.0, 0.0),
    N = 101,
    second_order = true,
    name = second_order ? "second order retraction" : "retraction",
    log_range = range(limits[1], limits[2]; length = N),
    plot = false,
    slope_tol = 0.1,
    error::Symbol = :none,
    window = nothing,
)
    Xn = X ./ norm(M, p, X) # normalize tangent direction
    # function for the directional derivative
    #
    T = exp10.(log_range)
    # points `p_i` to evaluate the error function at
    points = [exp(M, p, Xn, t) for t in T]
    approx_points = [retract(M, p, Xn, t, retraction_method) for t in T]
    errors = [distance(M, p, q) for (p, q) in zip(points, approx_points)]
    return prepare_check_result(
        log_range,
        errors,
        second_order ? 3.0 : 2.0;
        exactness_tol = exactness_tol,
        io = io,
        name = name,
        plot = plot,
        slope_tol = slope_tol,
        error = error,
        window = window,
    )
end

function plot_slope end

"""
    plot_slope(x, y;
        slope=2,
        line_base=0,
        a=0,
        b=2.0,
        i=1,
        j=length(x)
    )

Plot the result from the verification functions on data `x,y` with two comparison lines

1) `line_base` + t`slope`  as the global slope(s) the plot could have
2) `a` + `b*t` on the interval [`x[i]`, `x[j]`] for some (best fitting) comparison slope

!!! note
    This function has to be implemented for a certain plotting package.
    loading [Plots.jl](https://docs.juliaplots.org/stable/) provides a default implementation.
"""
plot_slope(x, y)

"""
    prepare_check_result(log_range, errors, slope;
        exactness_to = 1e3*eps(eltype(errors)),
        io = nothing
        name = "estimated slope",
        plot = false,
        slope_tol = 0.1,
        error = :none,
    )

Given a range of values `log_range`, with computed `errors`,
verify whether this yields a slope of `slope` in log-scale

Note that if the errors are below the given tolerance and the method is exact,
no plot is be generated,

# Keyword arguments

* `exactness_tol`: is all errors are below this tolerance, the verification is considered to be exact
* `io`:            provide an `IO` to print the result to
* `name`:          name to display in the plot title
* `plot`:          whether to plot the result, see [`plot_slope`](@ref)
  The plot is in log-log-scale. This is returned and can then also be saved.
* `slope_tol`:     tolerance for the slope (global) of the approximation
* `error`:         specify how to handle errors, `:none`, `:info`, `:warn`, `:error`
"""

function prepare_check_result(
    log_range,
    errors,
    slope;
    io::Union{IO,Nothing} = nothing,
    name = "estimated slope",
    slope_tol = 1e-1,
    plot = false,
    error::Symbol = :none,
    window = nothing,
    exactness_tol = 1e3 * eps(eltype(errors)),
)
    if max(errors...) < exactness_tol
        (io !== nothing) && print(
            io,
            "All errors are below the exactness tolerance $(exactness_tol). Your check can be considered exact, hence there is no use to check for a slope.\n",
        )
        return true
    end
    x = log_range[errors .> 0]
    T = exp10.(x)
    y = log10.(errors[errors .> 0])
    (a, b) = find_best_slope_window(x, y, length(x))[1:2]
    if isapprox(b, slope; atol = slope_tol)
        plot && return plot_slope(
            T,
            errors[errors .> 0];
            slope = slope,
            line_base = errors[1],
            a = a,
            b = b,
            i = 1,
            j = length(y),
        )
        (io !== nothing) && print(
            io,
            "Your $name's slope is globally $(@sprintf("%.4f", b)), so within $slope ± $(slope_tol).\n",
        )
        return true
    end
    # otherwise
    # find best contiguous window of length w
    (ab, bb, ib, jb) = find_best_slope_window(x, y, window; slope_tol = slope_tol)
    msg = "The $(name) fits best on [$(T[ib]),$(T[jb])] with slope  $(@sprintf("%.4f", bb)), but globally your slope $(@sprintf("%.4f", b)) is outside of the tolerance $slope ± $(slope_tol).\n"
    (io !== nothing) && print(io, msg)
    plot && return plot_slope(
        T,
        errors[errors .> 0];
        slope = slope,
        line_base = errors[1],
        a = ab,
        b = bb,
        i = ib,
        j = jb,
    )
    (error === :info) && @info msg
    (error === :warn) && @warn msg
    (error === :error) && throw(ErrorException(msg))
    return false
end

"""
    (a,b,i,j) = find_best_slope_window(X,Y,window=nothing; slope=2.0, slope_tol=0.1)

Check data X,Y for the largest contiguous interval (window) with a regression line fitting “best”.
Among all intervals with a slope within `slope_tol` to `slope` the longest one is taken.
If no such interval exists, the one with the slope closest to `slope` is taken.

If the window is set to `nothing` (default), all window sizes `2,...,length(X)` are checked.
You can also specify a window size or an array of window sizes.

For each window size, all its translates in the data is checked.
For all these (shifted) windows the regression line is computed (with `a,b` in `a + t*b`)
and the best line is computed.

From the best line the following data is returned

* `a`, `b` specifying the regression line `a + t*b`
* `i`, `j` determining the window, i.e the regression line stems from data `X[i], ..., X[j]`
"""
function find_best_slope_window(X, Y, window = nothing; slope = 2.0, slope_tol = 0.1)
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
