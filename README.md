# RandomBatchEwald

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://shiyibai5315.github.io/RBE.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://shiyibai5315.github.io/RBE.jl/dev/)
[![Build Status](https://github.com/shiyibai5315/RBE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/shiyibai5315/RBE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/shiyibai5315/RBE.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/shiyibai5315/RBE.jl)

### Introduction to the Random Batch Ewald (RBE) Method

The Random Batch Ewald (RBE) method is a novel algorithm designed for efficient molecular dynamics simulations of particle systems with long-range Coulomb interactions. Traditional methods for computing these interactions, such as the classical Ewald summation, often suffer from high computational complexity, typically \( \mathcal{O}(N^{3/2}) \) or \( \mathcal{O}(N \log N) \) per iteration, where \(N\) is the number of particles. 

In the RBE method, a random "minibatch" of Fourier modes is sampled at each simulation step, significantly reducing the computational burden while maintaining accuracy. This importance sampling reduces the variance of the force calculations, making the method both unbiased and efficient. The resulting algorithm achieves \( \mathcal{O}(N) \) complexity per time step, enabling simulations of large particle systems with long-range interactions that were previously computationally prohibitive.


### Mathematical Process of the RBE Method

The Random Batch Ewald (RBE) method is designed to efficiently calculate long-range Coulomb interactions in particle systems by utilizing a stochastic approach to sampling the Fourier space (k-space). The method involves the following steps:

1. **Ewald Splitting**: The Coulomb potential is split into a short-range real-space part and a long-range reciprocal-space part.
2. **Generating the k-space Vectors**: Create a set of reciprocal lattice vectors (k-vectors) within a specified cutoff.
3. **Computing Probabilities**: Calculate the probabilities associated with each k-vector based on a Gaussian distribution.
4. **Sampling from the k-space**: Use importance sampling to randomly select a subset of k-vectors.
5. **Force Calculation**: Compute the forces in both real space and reciprocal space using the sampled k-vectors.

#### Mathematical Derivation

**Ewald Splitting**:
The Coulomb potential $\frac{1}{r}$ is split into two parts:
```math
 \frac{1}{r} = \frac{\text{erfc}(\sqrt{\alpha} r)}{r} + \frac{\text{erf}(\sqrt{\alpha} r)}{r} 
 ```
where \(\alpha\) is a screening parameter.

**Generating k-space Vectors**:
The k-space vectors are generated such that their magnitude is within a cutoff \(k_c\):
```math
\mathbf{k} = \left( \frac{2\pi m_1}{L_x}, \frac{2\pi m_2}{L_y}, \frac{2\pi m_3}{L_z} \right)
```
where \( \|\mathbf{k}\| \leq k_c \).

**Computing Probabilities**:
For each k-vector \(\mathbf{k}\), the probability is given by:
```math
P(\mathbf{k}) = \frac{\exp\left(-\frac{\|\mathbf{k}\|^2}{4\alpha}\right)}{S} 
```
where \(S\) is a normalization constant.

**Sampling from k-space**:
Randomly sample \(p\) k-vectors based on their probabilities to form a minibatch.

**Force Calculation**:
The long-range force on particle \(i\) due to the k-space component is:
```math
\mathbf{F}_i^{\text{long}} = -\sum_{j=1}^{p} \frac{4\pi \mathbf{k}_j q_i}{V \|\mathbf{k}_j\|^2} \Im \left( \exp(-i \mathbf{k}_j \cdot \mathbf{r}_i) \rho_k[j] \right)
```

## Overview of the workflow in RBE process. 

#### Core Workflow of the RBE Method

1. **Initialization and Sampling**:
   - It uses importance sampling to randomly select a subset of k-vectors from the entire k-space. This subset will be used to approximate the long-range interactions.

2. **Update Charge Density Fourier Coefficients (\(\rho_k\))**:
   - The function updates the Fourier coefficients of the charge density based on the current positions of the particles and the sampled k-space vectors.

3. **Long-Range Force Calculation**:
   - Using the sampled k-space vectors and the updated \(\rho_k\), the function calculates the long-range forces acting on each particle. 

4. **Short-Range Force Calculation**:
   - The function identifies neighboring particles within a specified cutoff distance and computes the short-range forces between them.

#### Julia Code Overview

Here's a simplified view of the `update_acceleration!` function:

```julia
function update_acceleration!(interaction::RBEInteractions{T}, neighborfinder::ExTinyMD.AbstractNeighborFinder, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number}
    # Initialization
    boundary = sys.boundary
    samples = Vector{Point{3,T}}(undef, interaction.p)
    map!(x -> Point(T.(x)...), samples, sample(interaction.k_set, weights(interaction.prob), interaction.p))
    
    # Update charge density Fourier coefficients
    update_rho_k!(interaction, sys.n_atoms, interaction.p, samples, sys, info)
    update_finder!(neighborfinder, info)

    # Long-range force calculation
    for i in 1:sys.n_atoms
        force_long = calculate_F_long(i, interaction.V, interaction.p, interaction.rho_k, info, samples, sys)
        info.particle_info[i].acceleration += force_long / sys.atoms[i].mass
    end

    # Short-range force calculation
    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, boundary, interaction.cutoff)
        if iszero(dist_sq)
            continue
        end
        force_vector = calculate_F_short(coord_1, coord_2, interaction.α, i, j, sys)
        force_short = Point(force_vector...)
        info.particle_info[i].acceleration += force_short / sys.atoms[i].mass
        info.particle_info[j].acceleration -= force_short / sys.atoms[j].mass
    end
end
```

#### Mathematical Justification

- **Initialization**: \( \mathcal{O}(N) \)
- **Sampling k-space Vectors**: \( \mathcal{O}(p) \)
- **Updating \(\rho_k\)**: \( \mathcal{O}(p) \)
- **Long-Range Force Calculation**: \( N \times p = \mathcal{O}(pN) \)
- **Short-Range Force Calculation**: \( \mathcal{O}(N) \) (assuming a constant number of neighbors)

Hence, the total complexity is:
\[ \mathcal{O}(N) + \mathcal{O}(p) + \mathcal{O}(pN) + \mathcal{O}(pN) + \mathcal{O}(N) = \mathcal{O}(pN) \]
Since p is a constant. The actual complexity is \( \mathcal{O}(N) \).

### Results

The results from the Random Batch Ewald (RBE) method for molecular dynamics simulations have been thoroughly analyzed and visualized. Below are the key findings:

1. **Execution Time vs Steps**:
     This linear trend indicates the efficiency and scalability of the algorithm for large-time-scale simulations.
    ![999 atoms](/docs/src/execution.png)
2. **Memory Allocations vs Steps**:
     The memory allocation graph shows a linear increase in memory usage with the number of steps.
  ![999 atoms](/docs/src/memory.png)
1. **Computing costs and the number of atoms analysis**:
    Computation Time: As the number of atoms increases, the time required for each simulation step increases exponentially. This is likely due to the increased complexity of force calculations and the need to manage a larger number of particle interactions.
    Memory Usage: The memory required for the simulation also grows exponentially with the number of atoms. This is because larger systems require more data storage for particle positions, velocities, forces, and other related information.
    ![999 atoms](/docs/src/output.png)


2. **Absolute force error analysis**:
     Comparing RBE results with EwaldSummation result, the absolute error graph indicates the variation in error. While there is some fluctuation, the majority of the errors stay within a reasonable range, confirming the method's consistency in maintaining accuracy.
    ![999 atoms](/docs/src/force_diff.png)

3. **Mean Absolute Error (MAE)**:
     Two graphs display the MAE over steps. The results show that the MAE remains low and stable throughout most of the steps, with occasional spikes. These spikes could be attributed to specific configurations or numerical instabilities.
    ![999 atoms](/docs/src/mae_1.png)
    ![999 atoms](/docs/src/mae_2.png)

4. **Root Mean Square Error (RMSE)**:
     The RMSE plots show similar patterns to the MAE, with low values throughout the simulation steps and a few significant spikes. This consistency further supports the RBE method's accuracy and reliability.
    ![999 atoms](/docs/src/rmse_1.png)
    ![999 atoms](/docs/src/rmse_2.png)
5. **Benchmark Result**:
These result contains 300 atoms with charge of one step indicates that the RBE method can efficiently handle single-step computations with relatively low memory usage and allocation counts. In this code, the main factor affecting the calculation is the size of p.
```julia
BenchmarkTools.Trial: 10000 samples with 1 evaluation.
 Range (min … max):  196.300 μs …   4.742 ms  ┊ GC (min … max): 0.00% … 81.01%
 Time  (median):     221.900 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   289.533 μs ± 224.206 μs  ┊ GC (mean ± σ):  1.53% ±  3.23%

  ▆█▇▅▂▁                              ▁▂▂▂▂▁▁▁                  ▂
  █████████▇▆▅▄▃▄▄▁▁▃▁▃▁▁▁▁▃▁▁▃▁▃▃▁▁▁▇██████████▆▇▆▆▄▄▄▃▁▄▃▃▅▆▆ █
  196 μs        Histogram: log(frequency) by time       1.08 ms <

 Memory estimate: 3.59 KiB, allocs estimate: 117.
 ```


#### Recommendations and Future Work

To address the observed limitations and enhance the performance of the RBE method, the following strategies are recommended:

1. **Adaptive Sampling**:
   - Implement adaptive sampling techniques for k-space vectors to dynamically adjust the number of samples based on the system size and required accuracy. This could help balance computational cost and precision, especially for very large simulations.

2. **Memory Management**:
   - Develop and integrate more efficient memory management schemes. Techniques such as memory pooling, efficient data serialization, and compression can help manage the increased memory demands of larger systems.

3. **Error Analysis and Correction**:
   - Investigate the sources of error spikes in the force calculations and implement corrective measures. This could involve refining the numerical methods used in the RBE algorithm or incorporating error correction techniques to improve overall accuracy.

4. **Benchmarking and Profiling**:
   - Regularly benchmark and profile the simulation code to identify areas for improvement. Continuous performance monitoring and profiling can help detect inefficiencies and guide targeted optimization efforts.

## Notes
Please refer to https://github.com/shiyibai5315/ExTinyMD.jl/tree/dev to obtain the lattest version.

## References
1. Liang, J., Xu, Z. and Zhao, Y., 2021. Random-batch list algorithm for short-range molecular dynamics simulations. The Journal of Chemical Physics, 155(4). https://doi.org/10.1063/5.0056515
