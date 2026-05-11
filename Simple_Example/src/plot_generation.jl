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

function generate_error_comparison_plot(
    N::Vector{Int64}, 
    T::theoretic_data,
    error_data::Vector{<:Vector{<:Vector{<:Real}}};
    y_lim_lower::Real = -Inf, 
    y_lim_upper::Real = Inf,
    title::String = ""
)::Figure

    I = collect(0:N[3])

    # evaluate best-case data
    best_case_error_data = hcat(error_data[1]...)'
    best_case_error_mean = map(x -> mean(x), eachcol(best_case_error_data))
    best_case_error_standard_error =
        2 .* map(x -> std(x), eachcol(best_case_error_data)) ./ sqrt(N[2])

    # evaluate adaptive data
    adaptive_error_data = hcat(error_data[2]...)'
    adaptive_error_mean = map(x -> mean(x), eachcol(adaptive_error_data))
    adaptive_error_standard_error =
        2 .* map(x -> std(x), eachcol(adaptive_error_data)) ./ sqrt(N[2])

    # determine highest and lowest plot points
    best_case_lowest = best_case_error_mean .- best_case_error_standard_error
    adaptive_lowest = adaptive_error_mean .- adaptive_error_standard_error

    best_case_highest = best_case_error_mean .+ best_case_error_standard_error
    adaptive_highest = adaptive_error_mean .+ adaptive_error_standard_error

    # set x-axis details
    x_ticks = collect(range(0, N[3], length = 6))
    x_tick_labels = string.(round.(Int, x_ticks))

    Λ = T.A .* (2 ./ (2 .+ I)) .^ T.r

    value_collection_min = [best_case_lowest..., adaptive_lowest..., Λ...]
    value_collection_max = [best_case_highest..., adaptive_highest..., Λ...]

    y_minimum = 0.99 * minimum(value_collection_min)
    y_maximum = 1.01 * maximum(value_collection_max)

    if y_lim_lower > -Inf
        y_minimum = y_lim_lower
    end

    if y_lim_upper < Inf
        y_maximum = y_lim_upper
    end

    # for log10 scale we cannot have negative values
    y_minimum = max(y_minimum, eps(Float64))

    figure = Figure()

    axis = Axis(
        figure[1, 1],
        title = title,
        xlabel = "Number of Sample Points / Iterations",
        ylabel = "Absolute Error",
        yscale = log10,
        xticks = (x_ticks, x_tick_labels),
        limits = ((0, N[3]), (y_minimum, y_maximum))
    )

    band!(
        axis,
        I,
        max.(best_case_lowest, eps(Float64)),
        best_case_highest,
        color = (:blue, 0.2),
        label = nothing
    )

    lines!(
        axis,
        I,
        best_case_error_mean,
        linewidth = 2,
        color = :blue,
        label = "best-case error"
    )

    band!(
        axis,
        I,
        max.(adaptive_lowest, eps(Float64)),
        adaptive_highest,
        color = (:red, 0.2),
        label = nothing
    )

    lines!(
        axis,
        I,
        adaptive_error_mean,
        linewidth = 2,
        linestyle = :dot,
        color = :red,
        label = "adaptive error"
    )

    lines!(
        axis,
        I,
        Λ,
        linewidth = 2,
        linestyle = :dash,
        color = :green,
        label = "theoretical bound"
    )

    axislegend(axis, position = :rt)

    return figure
end

#----------------------------------------------------------------------------------#

function generate_accelerated_error_plot(
    N::Vector{Int64}, 
    T::theoretic_data,
    errors::Union{Nothing, Vector{<:Vector{<:Real}}};
    spacing::Int64 = 1, 
    y_lim_lower::Real = -Inf, 
    y_lim_upper::Real = Inf,
    title::String = ""
)::Figure

     I = collect(0:spacing:N[3])

    if errors === nothing
        adaptive_error_mean = Inf * ones(length(I))
        adaptive_error_standard_error = zeros(length(I))
    else
        # evaluate adaptive data
        adaptive_error_data = hcat(errors...)'
        adaptive_error_data = adaptive_error_data[:, I .+ 1]

        adaptive_error_mean = map(x -> mean(x), eachcol(adaptive_error_data))
        adaptive_error_standard_error =
            2 .* map(x -> std(x), eachcol(adaptive_error_data)) ./ sqrt(N[2])
    end

    # determine highest and lowest plot points
    adaptive_lowest = adaptive_error_mean .- adaptive_error_standard_error
    adaptive_highest = adaptive_error_mean .+ adaptive_error_standard_error

    # set x-axis details
    x_ticks = collect(range(0, N[3], length = 6))
    x_tick_labels = string.(round.(Int, x_ticks))

    Λ = T.A_plus .* (2 ./ (2 .+ I)) .^ T.r
    Λ_quadratic = T.B .* (2 ./ (1 .+ I)) .^ (2 * T.r)

    N_intersect = findfirst(Λ .> Λ_quadratic)

    if N_intersect === nothing
        Λ_optimal_combined = Λ
        Λ_linear_continued = nothing
    else
        Λ_optimal_combined = [Λ[1:(N_intersect - 1)]..., 
            Λ_quadratic[N_intersect:end]...]

        Λ_linear_continued = Λ[N_intersect:end]
    end

    value_collection_min = [adaptive_lowest..., Λ_optimal_combined...]
    value_collection_max = [adaptive_highest..., Λ_optimal_combined...]

    y_minimum = 0.99 * minimum(value_collection_min)
    y_maximum = 1.01 * maximum(value_collection_max)

    if y_lim_lower > -Inf
        y_minimum = y_lim_lower
    end

    if y_lim_upper < Inf
        y_maximum = y_lim_upper
    end

    # for log10 scale we cannot have negative values
    y_minimum = max(y_minimum, eps(Float64))

    figure = Figure(size = (1200, 800))

    axis = Axis(
        figure[1, 1],
        title = title,
        xlabel = "Number of Sample Points / Iterations",
        ylabel = "Absolute Error",
        yscale = log10,
        xticks = (x_ticks, x_tick_labels),
        limits = ((0, N[3]), (y_minimum, y_maximum))
    )

    band!(
        axis,
        I,
        max.(adaptive_lowest, eps(Float64)),
        adaptive_highest,
        color = (:red, 0.2),
        label = nothing
    )

    if errors === nothing
        adaptive_label = nothing
    else
        adaptive_label = "adaptive error"
    end

    lines!(
        axis,
        I,
        adaptive_error_mean,
        linewidth = 2,
        color = :red,
        label = adaptive_label
    )

    if N_intersect !== nothing
        lines!(
            axis,
            I[N_intersect:end],
            Λ_linear_continued,
            linewidth = 2,
            linestyle = :dot,
            color = :green,
            label = "default theoretical bound"
        )
    end

    lines!(
        axis,
        I,
        Λ_optimal_combined,
        linewidth = 2,
        linestyle = :dash,
        color = :green,
        label = "improved theoretical bound"
    )

    axislegend(axis, position = :rt)

    return figure
end

#----------------------------------------------------------------------------------#