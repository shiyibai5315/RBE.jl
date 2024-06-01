struct RBEInteraction{T} <: ExTinyMD.AbstractInteraction
    α::T
    cutoff::T
    ε::T
end

Base.show(io::IO, interaction::RBEInteraction) = print(io, "RBEInteraction with α = $(interaction.α), cutoff = $(interaction.cutoff), ε = $(interaction.ε)")

RBEInteraction(;α::T = 1.0, cutoff::T = 3.5, ε::T = 1.0) where T = RBEInteraction(α, cutoff, ε)



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


function generate_k_vector(α::Float64, L::Float64) #sampling for 3 times
    kx = mh_sample(α, L)
    ky = mh_sample(α, L)
    kz = mh_sample(α, L)
    return [kx, ky, kz]
end


function calculate_Fi(i::Int, p::Int, L::Float64, α::Float64, charges::Vector{Float64}, positions::Matrix{Tuple{Float64, Float64, Float64}})
    V = L^3
    S = calculate_S(α, L)
    Fi = zeros(Float64, 3)
    qi = charges[i]
    ri = [positions[i]...]
    
    kl = generate_k_vector(α,L) #sampling 3 times to have kl
    k2_ell = sum(kl .^ 2) #magnitude of kl
    rho_k = sum(charges[j] * exp(1im * dot(kl, [positions[j]...])) for j in 1:length(charges))

    for i in 1 : p
        exp_term = exp(-1im * dot(kl, ri)) #item in the force equation
        Fi += - (S / p) * ((4 * π * kl * qi) / (V * k2_ell)) * imag(exp_term * rho_k) #calculate force
    end
    
    return Fi
end


function ExTinyMD.energy()

end

function ExTinyMD.update_acceleration!(interaction::RBEInteraction, neighborfinder, sys, info)
    n_atoms = 2
    L = 50.0
    α = 3.0

    charges = [sys.atoms[i].charge for i in 1:n_atoms]
    positions = hcat([info.particle_info[i].position.coo for i in 1:n_atoms]...)
    update_finder!(neighborfinder, info)
    # calculate Fi
    # update the particle info according to Fi
    for i in 1:n_atoms
        Fi = RBE.calculate_Fi(i, n_atoms, L, α, charges, positions)
        info.particle_info[i].acceleration += Point{3, Float64}(Tuple(Fi / sys.atoms[i].mass))
    end
    
    return nothing
end
