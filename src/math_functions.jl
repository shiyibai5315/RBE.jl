function calculate_H(α::Float64, L::Float64; truncate::Bool=true, terms::Int=1)
    if truncate
        return sqrt(α * L^2 / π) * (1 + 2 * exp(-α * L^2))
    else
        H = 0.0
        for m in -terms:terms
            H += exp(-π^2 * m^2 / (α * L^2))
        end
        return sqrt(α * L^2 / π) * H
    end
end

# 计算S
function calculate_S(α::Float64, L::Float64)
    H = calculate_H(α, L)
    return H^3 - 1
end

# 计算概率
function calculate_probability(k::Vector{Float64}, α::Float64, S::Float64)
    k2 = sum(k .^ 2)
    return exp(-k2 / (4 * α)) / S
end

end # module MathFunctions