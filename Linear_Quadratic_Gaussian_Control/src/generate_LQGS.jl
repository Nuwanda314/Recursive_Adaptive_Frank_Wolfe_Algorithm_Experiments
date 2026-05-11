#----------------------------------------------------------------------------------#

function generate_extended_LQGS(
    d::Int64, 
    T::Int64
    )::LQGS_extended

    Id = diagm(ones(d))

    p = ceil(Int64, d / 2)
    Ip = diagm(ones(p))

    # generate A to have 1 one the diagonal and superdiagonal
    A = copy(Id)
    for i = 1:(d - 1)
        A[i, i + 1] = 1
    end
    A = [A for _ = 1:T]

    B = [Id for _ = 1:T]
    R = [Id / d ^ 2 for _ = 1:T]
    Q = [Id / d ^ 2 for _ = 1:(T + 1)]

    # generate C to project onto the first p components
    C = hcat(Ip, zeros(p, d - p))
    C = [C for _ = 1:T]

    X = Id
    W = [0.05 * Id for _ = 1:T]
    V = [0.1 * Ip for _ = 1:T]

    return extend_LQGS(LQGS(A, B, C, Q, R, X, W, V, T))
end

#----------------------------------------------------------------------------------#
S = generate_extended_LQGS(5, 5)