using RBE
using ExTinyMD
using Random

begin
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
    # #--------------------------MD process--------------------------
    n_atoms = 99
    atoms = create_atoms([(n_atoms ÷ 3, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 3, Atom(type = 2, mass = 1.0, charge = - 1.0)), (n_atoms ÷ 3, Atom(type = 3, mass = 10.0, charge = 0.0))])

    min_r = 0.1
    temp = 1.0
    L = 100.0
    boundary = CubicBoundary(100.0)
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary, min_r=min_r, max_attempts=100, rng=Random.GLOBAL_RNG, temp=temp)
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (RBEInteraction(), CellList3D(info, 4.5, boundary, 100))
    ]
    loggers = [TrajectoryLogger(;step = 10, trajectory_file="trajectory.txt",output = true), TemperatureLogger(100)]

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
    
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    simulate!(simulator, sys, info, 100)
    
end