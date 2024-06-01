using RBE
using ExTinyMD


#begin
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
    # charges = [1.0, -1.0, 0.5, -0.5, 1.5]
    # positions = [1.0 2.0 3.0 4.0 5.0; 0.0 -1.0 1.0 0.5 -0.5; 2.5 3.5 1.5 -2.0 0.0]
    
    i = 1  # Index of the atom for which we want to calculate the force

    # force = calculate_Fi(i, p, L, α, charges, positions)
    # println("Force on atom $i: $force")
    
    # #--------------------------MD process--------------------------
    n_atoms = 2
    atoms = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0, charge = 1.0))])

    min_r = 0.1
    temp = 1.0
    L = 50.0
    boundary = CubicBoundary(50.0)
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary, min_r=min_r, max_attempts=100, rng=Random.GLOBAL_RNG, temp=temp)
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (RBEInteraction(), CellList3D(info, 4.5, boundary, 100))
    ]
    loggers = Vector{ExTinyMD.AbstractLogger}()

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    # for (interaction, neighborfinder) in sys.interactions
    #     if isa(interaction, LennardJones)
    #         ExTinyMD.update_acceleration!(interaction, neighborfinder, sys, info)
    #     end
    # end

    # 计算额外的力 Fi 并添加到加速度中
    charges = [sys.atoms[i].charge for i in 1:n_atoms]
    positions = hcat([info.particle_info[i].position.coo for i in 1:n_atoms]...)


    for i in 1:n_atoms
        Fi = RBE.calculate_Fi(i, n_atoms, L, α, charges, positions)
        @show(Fi)
        info.particle_info[i].acceleration += Point{3, Float64}(Tuple(Fi / sys.atoms[i].mass))
    end
    

    println("Accelerations: ", [info.particle_info[i].acceleration for i in 1:n_atoms])
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    simulate!(simulator, sys, info, 10000)
    
#end







