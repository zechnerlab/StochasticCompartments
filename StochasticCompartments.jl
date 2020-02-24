module StochasticCompartments

using DifferentialEquations
using Distributions
using PyPlot
using StatsBase
using Random

include("SC_structs.jl")
include("SC_models.jl")
include("SC_SSA.jl")
include("SC_momentODEs.jl")
include("SC_figures.jl")

export System, TransitionClass
export add_transition_class, new_chemical_reaction_class
export SSA, solve_moment_equations
export figure_1, figure_birthdeath, figure_SI1
export figure_coagulationfragmentation
export figure_cellcommunication, figure_cellactivation
export figure_stemcellstart, figure_stemcelldynamics, figure_stemcellparameters, figure_stemcellperturbation

end  # END MODULE
