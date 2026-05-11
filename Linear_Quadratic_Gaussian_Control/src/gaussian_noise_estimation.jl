#----------------------------------------------------------------------------------#

function initial_estimations(
    system::LQGS_extended,
    N::Int64
    )::Vector{<:Matrix{<:Real}}

    # generate new sample points
    x = rand(MvNormal(zeros(system.n), system.X), N)'
    w = [rand(MvNormal(zeros(system.n), system.W[t]), N)' for t = 1:system.T]
    v = [rand(MvNormal(zeros(system.p), system.V[t]), N)' for t = 1:system.T]

    # compute covariance estimations
    X_estimated = (x' * x) / N
    W_estimated = [(w[t]' * w[t]) / N for t = 1:system.T]
    V_estimated = [(v[t]' * v[t]) / N for t = 1:system.T]

    return [X_estimated, W_estimated..., V_estimated...]
end

#----------------------------------------------------------------------------------#

function update_estimations(
    system::LQGS_extended, 
    Z_estimated::Vector{<:Matrix{<:Real}}, 
    N::Int64,
    k::Int64
    )::Vector{<:Matrix{<:Real}}

    # generate new sample points
    x = rand(MvNormal(zeros(system.n), system.X))
    w = [rand(MvNormal(zeros(system.n), system.W[t])) for t = 1:system.T]
    v = [rand(MvNormal(zeros(system.p), system.V[t])) for t = 1:system.T]

    α = 1 / (N + k)

    # recursively update covariance estimations
    Z_estimated[1] = (1 - α) * Z_estimated[1] + α * x * x'
    Z_estimated[2:(system.T + 1)] = (1 - α) * Z_estimated[2:(system.T + 1)] + 
        α * [w[t] * w[t]' for t = 1:system.T]
    Z_estimated[(system.T + 2):end] = (1 - α) * Z_estimated[(system.T + 2):end] + 
        α * [v[t] * v[t]' for t = 1:system.T]

    # ensure numerical symmetrie of each component of Z_estimated
    return [(Z_estimated[t] + Z_estimated[t]') / 2 for t = 1:(2 * system.T + 1)]
end

#----------------------------------------------------------------------------------#