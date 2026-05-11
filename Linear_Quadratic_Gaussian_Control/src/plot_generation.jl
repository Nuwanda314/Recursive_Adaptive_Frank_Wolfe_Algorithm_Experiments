#----------------------------------------------------------------------------------#

# set some basic plot properties beforehand
set_theme!(
    size = (1200, 800),
    figure_padding = (0, 40, 0, 0),

    Axis = (
        titlesize = 30,
        titlegap = 10,
        xlabelsize = 25,
        xlabelpadding = 15,
        ylabelsize = 25,
        ylabelpadding = 5,
        xticklabelsize = 15,
        yticklabelsize = 15,
        xgridvisible = true,
        ygridvisible = true,
    ),

    Legend = (
        patchsize = (30, 0),
        labelsize = 20,
        rowgap = 3,
        padding = (10, 10, 8, 8),
    ),

    Lines = (
        linewidth = 2,
    ),
)

#----------------------------------------------------------------------------------#

function generate_trajectory_comparison_plot(
    optimal_data::Vector{<:Real},
    trajectory_data::Vector{<:Vector{<:Vector{<:Real}}},
    N::Vector{Int64};
    y_lim_lower::Real = -Inf,
    y_lim_upper::Real = Inf,
)::Figure

    # evaluate best-case data
    best_case_trajectory_data = hcat(trajectory_data[1]...)'
    best_case_trajectory_mean = vec(map(mean, eachcol(best_case_trajectory_data)))
    best_case_trajectory_standard_error =
        vec(map(std, eachcol(best_case_trajectory_data)) ./ sqrt(N[2]))

    # evaluate adaptive data
    adaptive_trajectory_data  = hcat(trajectory_data[2]...)'
    adaptive_trajectory_mean  = vec(map(mean, eachcol(adaptive_trajectory_data)))
    adaptive_trajectory_standard_error =
        vec(map(std, eachcol(adaptive_trajectory_data)) ./ sqrt(N[2]))

    # determin highest and lowest plot points
    best_case_lowest = best_case_trajectory_mean .- best_case_trajectory_standard_error
    adaptive_lowest  = adaptive_trajectory_mean .- adaptive_trajectory_standard_error

    best_case_highest = best_case_trajectory_mean .+ best_case_trajectory_standard_error
    adaptive_highest  = adaptive_trajectory_mean .+ adaptive_trajectory_standard_error

    y_minimum = 0.99 * minimum([best_case_lowest; adaptive_lowest])
    y_maximum = 1.01 * maximum([best_case_highest; adaptive_highest])

    if y_lim_lower > -Inf
        y_minimum = y_lim_lower
    end

    if y_lim_upper < Inf
        y_maximum = y_lim_upper
    end

    x = 1:N[3]

    figure = Figure()

    axis = Axis(
        figure[1, 1],
        title = "(a) Trajectory Comparison",
        xlabel = "Number of Sample Points / Iterations",
        ylabel = "Cost of Control",
        limits = ((1, N[3]), (y_minimum, y_maximum)),
    )

    lines!(
        axis,
        x,
        optimal_data[1] .* ones(N[3]),
        color = :green,
        linestyle = :dashdot,
        label = "optimal value",
    )

    lines!(
        axis,
        x,
        optimal_data[2] .* ones(N[3]),
        color = :green,
        linestyle = :dash,
        label = "robust optimal value",
    )

    band!(
        axis,
        x,
        best_case_trajectory_mean .- best_case_trajectory_standard_error,
        best_case_trajectory_mean .+ best_case_trajectory_standard_error,
        color = (:blue, 0.25),
    )

    lines!(
        axis,
        x,
        best_case_trajectory_mean,
        color = :blue,
        label = "best-case trajectory",
    )

    band!(
        axis,
        x,
        adaptive_trajectory_mean .- adaptive_trajectory_standard_error,
        adaptive_trajectory_mean .+ adaptive_trajectory_standard_error,
        color = (:red, 0.25),
    )

    lines!(
        axis,
        x,
        adaptive_trajectory_mean,
        color = :red,
        linestyle = :dot,
        label = "adaptive trajectory",
    )

    axislegend(axis, position = :rc)

    return figure
end

#----------------------------------------------------------------------------------#

function generate_error_comparison_plot(
    error_data::Vector{<:Vector{<:Vector{<:Real}}},
    N::Vector{Int64};
    y_lim_lower::Real = -Inf,
    y_lim_upper::Real = Inf,
)::Figure

    # evaluate best-case data
    best_case_error_data = hcat(error_data[1]...)'
    best_case_error_mean = vec(map(mean, eachcol(best_case_error_data)))
    best_case_error_standard_error =
        vec(map(std, eachcol(best_case_error_data)) ./ sqrt(N[2]))

    # evaluate adaptive data
    adaptive_error_data  = hcat(error_data[2]...)'
    adaptive_error_mean  = vec(map(mean, eachcol(adaptive_error_data)))
    adaptive_error_standard_error =
        vec(map(std, eachcol(adaptive_error_data)) ./ sqrt(N[2]))

    # determin highest and lowest plot points
    best_case_lowest = best_case_error_mean .- best_case_error_standard_error
    adaptive_lowest  = adaptive_error_mean .- adaptive_error_standard_error

    best_case_highest = best_case_error_mean .+ best_case_error_standard_error
    adaptive_highest  = adaptive_error_mean .+ adaptive_error_standard_error

    y_minimum = 0.99 * minimum([best_case_lowest; adaptive_lowest])
    y_maximum = 1.01 * maximum([best_case_highest; adaptive_highest])

    if y_lim_lower > -Inf
        y_minimum = y_lim_lower
    end

    # for log10 scale we cannot have negative values
    y_minimum = max(y_minimum, eps(Float64))

    if y_lim_upper < Inf
        y_maximum = y_lim_upper
    end

    x = 1:N[3]

    figure = Figure()

    axis = Axis(
        figure[1, 1],
        title = "(b) Error Comparison",
        xlabel = "Number of Sample Points / Iterations",
        ylabel = "Absolute Error",
        limits = ((1, N[3]), (y_minimum, y_maximum)),
        yscale = log10,
    )

    band!(
        axis,
        x,
        max.(best_case_error_mean .- best_case_error_standard_error, eps(Float64)),
        best_case_error_mean .+ best_case_error_standard_error,
        color = (:blue, 0.2),
    )

    lines!(
        axis,
        x,
        best_case_error_mean,
        color = :blue,
        label = "best-case error",
    )

    band!(
        axis,
        x,
        max.(adaptive_error_mean .- adaptive_error_standard_error, eps(Float64)),
        adaptive_error_mean .+ adaptive_error_standard_error,
        color = (:red, 0.2),
    )

    lines!(
        axis,
        x,
        adaptive_error_mean,
        color = :red,
        linestyle = :dot,
        label = "adaptive error",
    )

    axislegend(axis, position = :rt)

    return figure
end

#----------------------------------------------------------------------------------#