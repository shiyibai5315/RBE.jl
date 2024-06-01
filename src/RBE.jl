module RBE

using SpecialFunctions
using Distributions
using LinearAlgebra
using ExTinyMD
using Random

include("math_functions.jl")
include("mh_sampler.jl")

export calculate_H, calculate_S, calculate_probability, mh_sample, calculate_force, calculate_Fi, generate_k_vector
export rbe_acceleration!


end
