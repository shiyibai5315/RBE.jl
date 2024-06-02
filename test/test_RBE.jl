using Test
using RBE
using Random
using LinearAlgebra

# Test for calculate_H
@testset "Test calculate_H" begin
    α = 3.0
    L = 50.0
    H = calculate_H(α, L)
    @test isapprox(H, 48.86, atol=0.01)  
end

# Test for calculate_S
@testset "Test calculate_S" begin
    α = 3.0
    L = 50.0
    S = calculate_S(α, L)
    @test isapprox(S, 116644.2574, atol=0.01)  
end

# Test for calculate_probability
@testset "Test calculate_probability" begin
    α = 3.0
    L = 50.0
    k = [1.0, 0.0, 0.0]
    S = calculate_S(α, L)
    probability = calculate_probability(k, α, S)
    @test isapprox(probability, 7.887610025789391e-6, atol=0.001)  
end

# Test for mh_sample
@testset "Test mh_sample" begin
    α = 3.0
    L = 50.0
    sample = mh_sample(α, L)
    @test sample != 0.0
end

# Test for generate_k_vector
@testset "Test generate_k_vector" begin
    α = 3.0
    L = 50.0
    k_vector = generate_k_vector(α, L)
    @test length(k_vector) == 3
end

# Test for calculate_Fi
@testset "Test calculate_Fi" begin
    α = 3.0
    L = 50.0
    p = 5
    charges = [1.0, -1.0]
    positions = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]
    positions_matrix = hcat([Tuple(position) for position in positions]...)
    rho_k = Complex{Float64}[0.1 + 0.1im, 0.2 + 0.2im, 0.3 + 0.3im, 0.4 + 0.4im, 0.5 + 0.5im]
    samples = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
    Fi = calculate_Fi(1, p, L, α, charges, positions_matrix, rho_k, samples)
    @test length(Fi) == 3
end

# Test for update_rho_k
@testset "Test update_rho_k" begin
    α = 3.0
    L = 50.0
    p = 5
    charges = [1.0, -1.0]
    positions = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]
    positions_matrix = hcat([Tuple(position) for position in positions]...)
    samples = Any[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
    rho_k = update_rho_k(p, charges, positions_matrix, samples)
    @test length(rho_k) == p
end

# Test for calculate_G
@testset "Test calculate_G" begin
    r = 1.0
    α = 3.0
    G = calculate_G(r, α)
    @test isapprox(G, 0.11161022509, atol=0.001)  
end

# Test for calculate_Fi_short
@testset "Test calculate_Fi_short" begin
    α = 3.0
    L = 50.0
    p = 5
    charges = [1.0, -1.0]
    positions = [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)]
    positions_matrix = hcat([Tuple(position) for position in positions]...)
    Fi_short = calculate_Fi_short(1, p, L, α, charges, positions_matrix)
    @test length(Fi_short) == 3
end