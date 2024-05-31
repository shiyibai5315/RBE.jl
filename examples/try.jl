using RBE



function example_usage()
    α = 3.0
    L = 50.0
    p = 5
    
    H = calculate_H(α, L)
    S = calculate_S(α, L)
    
    println("H = $H, S = $S")
    

    samples = mh_sample(α, L)
    println("Samples: ", samples)
    
    
    k = [1.0, 0.0, 0.0]
    probability = calculate_probability(k, α, S)
    println("Probability for k = $k: $probability")
 
    charges = [1.0, -1.0, 0.5, -0.5, 1.5]
    positions = [1.0 2.0 3.0 4.0 5.0; 0.0 -1.0 1.0 0.5 -0.5; 2.5 3.5 1.5 -2.0 0.0]
    
    i = 1  # Index of the atom for which we want to calculate the force

    force = calculate_Fi(i, p, L, α, charges, positions)
    println("Force on atom $i: $force")
    
    
end

example_usage()






