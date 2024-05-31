function calculate_H(α, L)
  
    const_part = sqrt(α * L^2 / π)
    
    exp_term(m, α, L) = exp(-α * m^2 * L^2)
    
    sum_terms = 0.0
    for m in -1:1
        sum_terms += exp_term(m, α, L)
    end
    
    H = const_part * sum_terms
    
    return H
end

function calculate_S(α::Float64, L::Float64)
    H = calculate_H(α, L)
    return H^3 - 1
end


function calculate_probability(k::Vector{Float64}, α::Float64, S::Float64)
    k2 = sum(k .^ 2)
    return exp(-k2 / (4 * α)) / S
end


function calculate_Fi(i::Int, p::Int, L::Float64, α::Float64, charges::Vector{Float64}, positions::Matrix{Float64})
    V = L^3
    S = calculate_S(α, L)
    Fi = zeros(Float64, 3)
    
    k_samples = [rand(Normal(0, sqrt(α * L^2 / (2 * π^2))), 3) for _ in 1:p] #produce the number of k samples
    
    qi = charges[i]
    ri = positions[:, i]
    
    for k_ell in k_samples
        k2_ell = sum(k_ell .^ 2)
        
        rho_k = sum(charges[j] * exp(1im * dot(k_ell, positions[:, j])) for j in 1:length(charges))
        exp_term = exp(-1im * dot(k_ell, ri)) #item in the summation equation

        Fi += - (S / size(k_samples, 1)) * ((4 * π * k_ell * qi) / (V * k2_ell)) * imag(exp_term * rho_k)
    end
    
    return Fi
end

