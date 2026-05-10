#----------------------------------------------------------------------------------#

struct LQGS
    # system evolution matrices
    A::Vector{<:Matrix{<:Real}} # state
    B::Vector{<:Matrix{<:Real}} # control
    C::Vector{<:Matrix{<:Real}} # observation

    # cost matrices
    Q::Vector{<:Matrix{<:Real}} # state
    R::Vector{<:Matrix{<:Real}} # control

    # gaussian noise matrices
    X::Matrix{<:Real}           # initial
    W::Vector{<:Matrix{<:Real}} # state
    V::Vector{<:Matrix{<:Real}} # observation

    # finite time horizon
    T::Int64
end

#----------------------------------------------------------------------------------#

struct LQGS_extended
    # system evolution matrices
    A::Vector{<:Matrix{<:Real}} # state
    B::Vector{<:Matrix{<:Real}} # control
    C::Vector{<:Matrix{<:Real}} # observation

    # cost matrices
    Q::Vector{<:Matrix{<:Real}} # state
    R::Vector{<:Matrix{<:Real}} # control

    # gaussian noise matrices
    X::Matrix{<:Real}           # initial
    W::Vector{<:Matrix{<:Real}} # state
    V::Vector{<:Matrix{<:Real}} # observation

    # finite time horizon
    T::Int64

    # include additional implicit information to minimize redundant computations and 
    # enhance performance:

    # space dimensions 
    n::Int64 # state
    m::Int64 # control
    p::Int64 # observation

    # roots of the semidefinite gaussian noise matrices
    root_X::Matrix{<:Real}           # initial
    root_W::Vector{<:Matrix{<:Real}} # state
    root_V::Vector{<:Matrix{<:Real}} # observation

    # discrete-time finite-horizon backwards Ricatti matrices
    P::Vector{<:Matrix{<:Real}}
end

#----------------------------------------------------------------------------------#

function extend_LQGS(
    system::LQGS
    )::LQGS_extended

    # extract space dimension out of the system evolution matrices
    n = size(system.A[1])[1] # state
    m = size(system.B[1])[2] # control
    p = size(system.C[1])[1] # observation

    # compute roots of the gaussian noise matrices
    root_X = √(system.X)                         # initial
    root_W = [√(system.W[t]) for t = 1:system.T] # state
    root_V = [√(system.V[t]) for t = 1:system.T] # observation

    P = Vector{Matrix{Float64}}(undef, system.T + 1)

    # compute finite-horizon discrete-time backwards Ricatti matrices
    P[system.T + 1] = system.Q[system.T + 1] # initialize backwards recursion
    for t = system.T:-1:1
        c₁ = system.B[t]' * P[t + 1]
        c₂ = c₁ * system.A[t]
        
        # compute inverse via \ operator
        c₃ = (system.R[t] + c₁ * system.B[t]) \ I
        
        P[t] = system.A[t]' * P[t + 1] * system.A[t] + system.Q[t] - c₂' * c₃ * c₂
    end

    return LQGS_extended(system.A, system.B, system.C, system.Q, system.R, system.X,
        system.W, system.V, system.T, n, m, p, root_X, root_W, root_V, P)
end

#----------------------------------------------------------------------------------#