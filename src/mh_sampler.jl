function acceptance_probability(m_star::Int, α::Float64, L::Float64)
    if m_star == 0
        return erf(1 / (2 * sqrt(α * L^2 / π^2)))
    else
        sqrt_val = sqrt(α * L^2 / π^2)
        return 0.5 * (erf((abs(m_star) + 0.5) / sqrt_val) - erf((abs(m_star) - 0.5) / sqrt_val))
    end
end

# MH samling
function mh_sample(α::Float64, L::Float64, num_samples::Int)
    samples = Int[]
    for _ in 1:num_samples
        x_star = rand(Normal(0, sqrt(α * L^2 / (2 * π^2))))
        m_star = round(Int, x_star)
        q = acceptance_probability(m_star, α, L)
        if rand() < q
            push!(samples, m_star)
        end
    end
    return samples
end
