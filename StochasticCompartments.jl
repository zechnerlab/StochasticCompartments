module StochasticCompartments

using DifferentialEquations
using Distributions
using PyPlot
using StatsBase
using Random
using SpecialFunctions
using HypergeometricFunctions

include("SC_structs.jl")
include("SC_models.jl")
include("SC_SSA.jl")
include("SC_momentODEs.jl")
include("SC_figures.jl")

export System, TransitionClass
export add_transition_class, new_chemical_reaction_class
export SSA, solve_moment_equations
export Figure_1
export Figure_2BC, Figure_2EF, Figure_S1
export Figure_3B, Figure_3C, Figure_3C_inset
export Figure_3E, Figure_3F, Figure_3G, Figure_3H, Figure_S4

end  # END MODULE
