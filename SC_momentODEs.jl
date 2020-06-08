
function solve_moment_equations(S::System,n0::Matrix{Int64},tspan::Tuple{Float64,Float64};Tpoints::Int64=91)
assert_model(S,n0)
@assert S.init_moments != nothing "Field init_moments is empty! Please provide a function in S.init_moments"
return solve_moment_equations(S,S.init_moments(n0),tspan;Tpoints=Tpoints)
end

function solve_moment_equations(S::System,Mom0::Vector{Float64},tspan::Tuple{Float64,Float64};Tpoints::Int64=91)
    tplot=collect(range(tspan[1],stop=tspan[2],length=Tpoints))
    Moments=zeros(length(Mom0),Tpoints)
    solve_moment_equations!(S,Mom0,tplot,Moments)
    return tplot,Moments
end


function solve_moment_equations!(S::System,Mom0::Vector{Float64},timepoints::Vector{Float64},Moments::Matrix{Float64})
    @assert S.moment_equations != nothing "Field moment_equations is empty! Please provide a function in S.moment_equations"
    prob=ODEProblem(S.moment_equations,Mom0,(timepoints[1],timepoints[end]),S)
    sol=solve(prob,alg_hints=[:stiff],saveat=timepoints)
    for t=1:length(timepoints)
        Moments[:,t] .= sol(timepoints[t]) # sol[t]
    end
end


# Moment Equations for "Nested Birth-Death Process"

function birth_death_ODEs(dM::Vector{Float64},M::Vector{Float64},S::System,t::Float64)
kb = S.c[1].k
kd = S.c[2].k
kB = S.c[3].k
lambda = S.c[3].parameters
kE = S.c[4].k
# N (Number of Compartments)
dM[1]= kB -kE*M[1]
# M (Total Mass)
dM[2]= kB*lambda -kE*M[2] +kb*M[1] -kd*M[2]
# S (Sum of squared content)
dM[3]= kB*lambda*(1+lambda) -kE*M[3] +kb*(M[1]+2*M[2]) +kd*(M[2]-2*M[3])
# N^2
dM[4]= kB*(1+2*M[1]) +kE*(M[1]-2*M[4])
#M^2
dM[5]= kB*(lambda*(1+lambda)+2*lambda*M[2]) +kE*(M[3]-2*M[5]) +kb*(M[1]+2*M[6]) +kd*(M[2]-2*M[5])
# N*M
dM[6]= kB*(lambda*(1+M[1])+M[2]) +kE*(M[2]-2*M[6]) +kb*M[4] -kd*M[6]
return
end

function birth_death_initial(n0::Matrix{Int64}) :: Vector{Float64}
M0=zeros(6)
M0[1] = size(n0,2)
M0[2] = sum(n0)
M0[3] = sum(n0.^2)
M0[4] = M0[1]^2
M0[5] = M0[2]^2
M0[6] = M0[1]*M0[2]
return M0
end


# Moment Equations for "Stochastic Coagulation-Fragmentation Dynamics"

function coagulation_fragmentation_ODEs(dM::Vector{Float64},M::Vector{Float64},S::System,t::Float64)
kC = S.c[1].k
kF = S.c[2].k
kI = S.c[3].k
lambda = S.c[3].parameters
kE = S.c[4].k
## Moments to close
XXX = 2*M[3]^2/M[2]-M[2]*M[3]/M[1] # Gamma for <X^3> ## 3*M[3]*M[2]/M[1] - 2*M[2]^3/M[1]^2 # M[1]*(M[3]/M[2])^3
N3 = 2*M[4]^2/M[1]-M[1]*M[4]  # Gamma for <N^3>
N2M = 2*M[4]*M[6]/M[1]-M[4]*M[2]  # Gamma for <N^2*M>
# N (Number of Compartments)
dM[1]= kI -kE*M[1] -kC/2*(M[4]-M[1]) +kF*M[2]
# M (Total Mass)
dM[2]= kI*lambda -kE*M[2]
# S (Sum of squared content)
dM[3]= kI*lambda*(1+lambda) -kE*M[3] +kC*(M[5]-M[3]) +kF/3*(M[3]-XXX)
# N^2
dM[4]= kI+2*kI*M[1] +kE*(M[1]-2*M[4]) +kC/2*(M[4]-M[1])-kC*(N3-M[4]) +kF*(M[2]+2*M[6])
# M^2
dM[5]= kI*lambda*(1+lambda+2*M[2]) +kE*(M[3]-2*M[5])
# N*M
dM[6]= kI*(lambda*(1+M[1])+M[2]) +kE*(M[2]-2*M[6]) +kC/2*(M[6]-N2M) +kF*M[5]
return
end

function coagulation_fragmentation_initial(n0::Matrix{Int64}) :: Vector{Float64}
M0=zeros(6)
M0[1] = size(n0,2)
M0[2] = sum(n0)
M0[3] = sum(n0.^2)
M0[4] = M0[1]^2
M0[5] = M0[2]^2
M0[6] = M0[1]*M0[2]
return M0
end


# Moment Equations for "Transciption dynamics in a cell community"

function cell_community_ODEs(dM::Vector{Float64},M::Vector{Float64},S::System,t::Float64)
kcom = S.c[1].k
kbG = S.c[2].k
kdG = S.c[3].k
kS = S.c[4].k
kbS = S.c[5].k
kdS = S.c[6].k
## Moments to close
MG3 = 2*M[4]^2/M[2]-M[2]*M[4]    #  Gamma <MG^3>
MG2MP = 2*M[4]*M[6]/M[2]-M[4]*M[3]  #  Gamma for <MG^2*MP>
# N (Number of Compartments)
dM[1]= 0
# MG (Total G Mass)
dM[2]= kcom*(M[1]*M[2]-M[4]) +kbG*(M[1]-M[2]) -kdG*M[2]
# MP (Total P Mass)
dM[3]= kbS*M[1] +kS*M[2] -kdS*M[3]
# MG^2
dM[4]= kcom*(M[1]*M[2]-M[4])+2*kcom*(M[1]*M[4]-MG3)
dM[4]+= kbG*(M[1]-M[2]+2*(M[1]*M[2]-M[4])) +kdG*(M[2]-2*M[4])
# MP^2
dM[5]= kbS*M[1]*(1+2*M[3]) +kS*(M[2]+2*M[6]) +kdS*(M[3]-2*M[5])
# MG*MP
dM[6]= kcom*(M[1]*M[6]-MG2MP) +kbS*M[1]*M[2] +kS*M[4]
dM[6]+= kbG*(M[1]*M[3]-M[6]) -(kdG+kdS)*M[6]
return
end

function cell_community_initial(n0::Matrix{Int64}) :: Vector{Float64}
M0=zeros(6)
M0[1] = size(n0,2)
M0[2] = sum(n0[1,:])
M0[3] = sum(n0[2,:])
M0[4] = M0[2]^2
M0[5] = M0[3]^2
M0[6] = M0[2]*M0[3]
return M0
end


# Moment Equations for "Stem Cell Population Dynamics"

function stemcells_ODEs(dM::Vector{Float64},M::Vector{Float64},S::System,t::Float64)
knf = S.c[1].k
kS = S.c[2].k
kF_asym = S.c[3].k
kF_sym  = S.c[4].k
kE = S.c[5].k
## Moments to close
GPPP = 2*M[4]^2/M[3]-M[3]*M[4]/M[2] # Gamma
NMG = M[1]*M[2]     # MF
NMGP = M[1]*M[3]    # MF
MGMGP = M[2]*M[3]   # MF
MGMGPP = M[2]*M[4]  # MF
MG3 = 2*M[6]^2/M[2]-M[2]*M[6]    # Gamma
# N (Number of Compartments)
dM[1]= (kF_asym+kF_sym)*M[3] -kE*(M[1]-M[2])
# M_G
dM[2]= kF_sym*M[3] -knf/2*(M[6]-M[2])
# M_GP
dM[3]= -(kF_asym+kF_sym)*M[4] +kS*M[2] -knf/2*(MGMGP-M[3])
# M_GPP
dM[4]= -(kF_asym+kF_sym)*GPPP +kS*(M[2]+2*M[3]) -knf/2*(MGMGPP-M[4])
# N^2
dM[5]=  (kF_asym+kF_sym)*(M[3]+2*NMGP) +kE*(M[1]-M[2]-2*(M[5]-NMG))
# M_G^2
dM[6]= kF_sym*(M[3]+2*MGMGP) +knf/2*(M[6]-M[2]) -knf*(MG3-M[6])
return
end

function stemcells_initial(n0::Matrix{Int64}) :: Vector{Float64}
M0=zeros(6)
M0[1] = size(n0,2)
M0[2] = sum(n0[1,:])
M0[3] = sum(n0[1,:] .* n0[2,:])
M0[4] = sum(n0[1,:] .* n0[2,:].^2)
M0[5] = M0[1]^2
M0[6] = M0[2]^2
return M0
end
