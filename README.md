# StochasticCompartments.jl

Julia code to perform stochastic simulations and solve moment equations for stochastic models of interacting biochemical compartments. 

The provided code reproduces the case studies shown in the paper "Stochastic reaction networks in dynamic compartment populations" (https://arxiv.org/abs/2002.00931) and can be also used to implement further models, upon minor coding efforts.

For further information, contact Lorenzo Duso (duso@mpi-cbg.de) & Christoph Zechner (zechner@mpi-cbg.de)


# Installation and requirements

Julia version greater or equal than 1.1 and the packages Distributions, PyPlot, DifferentialEquations, StatsBase, Random, SpecialFunctions and HypergeometricFunctions are required. 


# Reproducing the figures of the paper

To get started, include the module StochasticCompartments.jl
```julia
include("./StochasticCompartments.jl")
using .StochasticCompartments
```
The following functions are dedicated to reproduce the corresponding figures of the paper: 
```julia
Figure_1()                         # for the panels of Fig. 1 (concept figure)
Figure_2BC()                       # for Figs. 2B and 2C (nested birth-death model)
Figure_S1()                        # for Fig. S.1
Figure_2EF()                       # for Figs. 2E, 2F and S.2 (coagulation-fragmentation case study)
Figure_3B()                        # for Figs. 3B and S.3 (dynamics of the protein expression)
Figure_3C()                        # for Figs. 3C and and S.3  (variance-to-mean ratio)
Figure_3C_inset()                  # for the inset of Fig. 3C (steady-state protein distributions)
Figure_3E()                        # for Fig. 3E (stem cell lineage realization)
Figure_S4()                        # for Fig. S4 (division time distribution)
Figure_3F()                        # for Fig. 3F (stem cell moment dynamics)
Figure_3G()                        # for Figs. 3G and S.5 (stem cell steady state)
Figure_3H()                        # for Fig. 3H (stem cell perturbation)
```
Note that some of figures show the output of single stochastic realizations, but the random number generator is reset always to the same condition to ensure perfect reproducibility. In order to obtain different random realizations, it is sufficient to pass the keyword `seed = nothing`, for instance `Figure_1(seed=nothing)` .

# Advanced usage

MODEL DECLARATION 

This section exaplains how to define new models and independently perform simulations.

The user can define a model through a variable of type `System`, which comprises the following fields:
```julia
    name::String
    n_species::Int64                                  # Number of chemical species
    transition_classes::Vector{TransitionClass}       # Defining model dynamics (used in SSA simulations)
    MomDict::Dict{Int64,Vector{Int64}}                # associates moment indices with exponents
    moment_equations::Union{Function,Nothing}         # Moment Equations f(dM,M,S,t)
    init_moments::Union{Function,Nothing}             # Based on implementation of Moment Equations
```
Note that the fields `moment_equations` and `init_moments` have to be implemented manually (for instance, as shown in the file SC_momentODEs.jl). The population dynamics is specified by the vector `c` of transition classes, whose implementation is related to that of `MomDict`. The latter is a dictionary that links some integer indices to some integer exponent vectors of length `n_species`. This allows the user to declare which moments the simulator needs to track to perform the simulation and/or that are desired in output. 
 
For instance, let's initialize a new model named `MyModel` comprising e.g. 2 chemical species. Let's assume that we are interested in following the dynamics of the total number of compartments (i.e. exponents `[0,0]`), the masses of species 1 and 2 (i.e. `[1,0]` and `[0,1]`) and their crossmoment `[1,1]`. We can assign these model settings to the variable 'S' as follows
```julia
S = System("MyModel",2, Dict(1=>[0,0],2=>[1,0],3=>[0,1],4=>[1,1]))
```
When designing the functions for each transition class associated with the model `S`, the indexing of the moments will follow the provided dictionary `MomDict`. Note, however, that the moment index `1` must always correspond to the total number of compartments.
Each transition class is stored in a variable of type `TransitionClass`, which comprises the following fields:
```julia
    rc::Int64                  # Number of reactant compartments
    pc::Int64                  # Number of product compartments
    DeltaN::Int64              # = pc - rc (just for conveniency)
    k::Float64                 # Rate constant
    parameters::Union{Nothing,Float64,Vector{Int64}}
    g::Union{Function,Nothing}         # Reactant compartments kernel
    pi::Union{Function,Nothing}        # Product compartments distribution
    H::Union{Function,Nothing}         # Total propensity
    fast_sample_reactants!::Union{Function,Nothing} # optional, optimized for specific class
```
A new transition class with `rc` reactant compartment, `pc` product compartments and rate constant `k` is initialized by running 
```julia
MyTransitionClass = TransitionClass(rc,pc,k)
```
Next, for each transition class it is necessary to declare:
1) the total class propensity function `H`, with arguments `(n::Matrix{Int64},Mom::Vector{Int64})`
2) the reactant kernel function `g`, with arguments `xc::Vector{Vector{Int64}}` OR the function `fast_sample_reactants!`. The latter can be manually implemented by the user and it serves to boost simulation speed (highly recommended). 
3) the outcome distribution, with arguments  `(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})` or `(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, param::Float64)`
where the integer matrix `n` has `n_species` rows and stores the state of the system by associating each compartment with one of its columns.

A customized function is provided to facilitate the declaration of transition classes associated with single-compartment chemical reactions. A chemical event with state-change vector `change_vector::Vector{Int64}` of length `n_species` and rate constant `k::Float64` is initialized by
```julia
MyChemicalClass = new_chemical_reaction_class(change_vector, k)
```
Afterwards, `MyChemicalClass` still requires the total class propensity `H` and the function `g` or `fast_sample_reactants!`.

To equip the variable `S` with some well-defined transition classes, it is sufficient to use the function 
```julia
add_transition_class(S,c...)
```
We invite the reader to examine the declaration of the case studies in the file `SC_models.jl` to see in practice how to specify models.

MODEL SIMULATION

The function `SSA` allows the user to perform stochastic simulations of a model `S` starting with initial condition `n0` of the compartment population. It returns different kinds of output depending on the type of arguments passed.

```julia
SSA(S::System, n0::Matrix{Int64}, timepoints::Vector{Float64})
```
performs one stochastic simulation from initial time `timepoints[1]` to final time `timepoints[end]`. It returns the final state of the system and the values of the population moments speficied in `S.MomDict` at each time in `timepoints`. If the keyword argument `full_story=true` is given, `SSA` returns the state of the system at each time in `timepoints` instead of only the final one.

```julia
SSA(S::System, n0::Matrix{Int64}, tmax::Float64)
```
performs one stochastic simulation from initial time `t=0` to final time `tmax`. It returns the final state of the system, the vector of event times occurred in the interval `[0,tmax)` and the population moments speficied in `S.MomDict` at each event. If the keyword argument `full_story=true` is given, it returns the state of the system at each event occurred in the interval `[0,tmax)` instead of only the final one.


```julia
SSA(S::System, n0::Matrix{Int64}, timepoints::Vector{T}, Nsamples::Int64) where T <: Real
```
performs `Nsamples` stochastic simulations from initial time `timepoints[1]` to final time `timepoints[end]`. It returns the Monte Carlo average and variance of the population moments speficied in `S.MomDict` at each time in `timepoints`. 

If the user has implemented the moment equations functions for the model and included them in the respective fields of the variable `S::System`, moment equations can be solved by the function
```julia
solve_moment_equations(S::System,n0::Matrix{Int64},tspan::Tuple{Float64,Float64} ; Tpoints::Int64)
```
where the tuple `tspan` specifies the initial and final simulation time. The keyword argument `Tpoints` can be used to assign at how many equally spaced points in the interval `tspan` the solution should be returned. The function returns a time vector and the corresponding value of the population moments.

Further clarifications can be found by reading the commented code or by contacting the authors at duso@mpi-cbg.de or zechner@mpi-cbg.de . 


