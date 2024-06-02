function acceptance_probability(m_star::Float64, α::Float64, L::Float64)
    if m_star == 0
        return erf(1 / (2 * sqrt(α * L^2 / π^2)))
    else
        sqrt_val = sqrt(α * L^2 / π^2)
        return 0.5 * (erf((abs(m_star) + 0.5) / sqrt_val) - erf((abs(m_star) - 0.5) / sqrt_val))
    end
end

function mh_sample(α::Float64, L::Float64)
    sample = 0.0
    while sample == 0.0
        x_star = rand(Normal(0, sqrt(α * L^2 / (2 * π^2))))
        m_star = Float64(round(Int, x_star))
        q = acceptance_probability(m_star, α, L)
        if rand() < q && m_star !=0
            sample =  m_star
            break
        end
    end
    return sample
end

function sampling(α, L, n)
    k_set = []
    for i in 0:2*L
        k_vector = generate_k_vector(α,L)
        push!(k_set, k_vector)
    end
    k_set = collect(k_set)
    z = sum(exp(-(norm(k)^2) / (4*α)) for k in k_set)
    prob = [exp(-(norm(k)^2) / (4*α)) / z for k in k_set]
    samples = sample(k_set, weights(prob), n)

    return samples
end