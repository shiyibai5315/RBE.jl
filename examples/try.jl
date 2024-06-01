using RBE
using ExTinyMD


function example_usage()
    α = 3.0
    L = 50.0
    p = 5
    
    #calculation
    H = calculate_H(α, L)
    S = calculate_S(α, L)
    
    println("H = $H, S = $S")
    

    samples = mh_sample(α, L)
    println("Samples: ", samples)
    
    
    k = [1.0, 0.0, 0.0]
    probability = calculate_probability(k, α, S)
    println("Probability for k = $k: $probability")
 
    #rbe process
    charges = [1.0, -1.0, 0.5, -0.5, 1.5]
    positions = [1.0 2.0 3.0 4.0 5.0; 0.0 -1.0 1.0 0.5 -0.5; 2.5 3.5 1.5 -2.0 0.0]
    
    i = 1  # Index of the atom for which we want to calculate the force

    force = calculate_Fi(i, p, L, α, charges, positions)
    println("Force on atom $i: $force")
    
    #MD process
    n_atoms = 5
    atoms = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0, charge = 1.0))])

    min_r = 0.1
    temp = 1.0
    L = 50.0
    boundary = CubicBoundary(50.0)
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary, min_r=min_r, max_attempts=100, rng=Random.GLOBAL_RNG, temp=temp)
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    interactions = Vector{Tuple{ExTinyMD.AbstractInteraction, ExTinyMD.AbstractNeighborFinder}}()
    loggers = Vector{ExTinyMD.AbstractLogger}()

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    rbe_acceleration!(sys, info, boundary; α)
    println("Accelerations: ", [info.particle_info[i].acceleration for i in 1:n_atoms])
end

example_usage()





