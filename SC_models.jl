# functions defining model of each case study

function ConceptExample()
    S=System("ConceptExample",2,Dict(1=>[0,0],2=>[1,0],3=>[0,1],4=>[2,0]))
    kI, lambda, kC, kr, kd = 1.0, 10.0, 0.01, 0.1, 0.1

    intake=TransitionClass(0,1,kI)
    intake.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> 1
    intake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, param::Float64)
                        yc[1][1]=rand(Poisson(param))
                        yc[1][2]=0
                end
    intake.parameters = lambda

    coag=TransitionClass(2,1,kC)
    coag.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]*(Mom[1]-1)/2
    coag.fast_sample_reactants! = fast_sample_uniform_2
    coag.pi = (yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}) -> for d=1:length(yc[1]) yc[1][d] = xc[1][d] + xc[2][d] end

    biconv=new_chemical_reaction_class([-2;1],kr)
    biconv.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> (Mom[4]-Mom[2])/2
    biconv.g = xc::Vector{Vector{Int64}} -> xc[1][1]*(xc[1][1]-1)/2

    death=new_chemical_reaction_class([0;-1],kd)
    death.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[3]
    death.g = xc::Vector{Vector{Int64}} -> xc[1][2]
    death.fast_sample_reactants! = fast_sample_mass_2

    add_transition_class(S,intake,coag,biconv,death)
    return S
end


# CASE STUDY "Nested Birth-Death Process"
function BirthDeath()
    S=System("BirthDeath",1,Dict(1=>[0],2=>[1]))
    kb, kd, kI, kE = 1.0, 0.1, 1.0, 0.01

    birth=new_chemical_reaction_class([1],kb)
    birth.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]
    birth.g = xc::Vector{Vector{Int64}} -> 1
    birth.fast_sample_reactants! = fast_sample_uniform_1

    death=new_chemical_reaction_class([-1],kd)
    death.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]
    death.g = xc::Vector{Vector{Int64}} -> xc[1][1]
    death.fast_sample_reactants! = fast_sample_mass_1

    intake=TransitionClass(0,1,kI)
    intake.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> 1
    intake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, param::Float64)
                    for d=1:length(yc[1])
                        yc[1][d]=rand(Poisson(param))
                    end
                end
    intake.parameters = birth.k / death.k

    exit=TransitionClass(1,0,kE)
    exit.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]
    exit.g = xc::Vector{Vector{Int64}} -> 1
    exit.fast_sample_reactants! = fast_sample_uniform_1

    add_transition_class(S,birth,death,intake,exit)

    S.moment_equations = birth_death_ODEs
    S.init_moments = birth_death_initial
    return S
end


# CASE STUDY "Stochastic Coagulation-Fragmentation Dynamics"

function CoagulationFragmentation()
    S=System("CoagulationFragmentation",1,Dict(1=>[0],2=>[1],3=>[2]))  # also Mom[2] just to look at it
    kC, kF, kI, kE = 0.005, 0.005, 10.0, 0.1

    coag=TransitionClass(2,1,kC)
    coag.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> div(Mom[1]*(Mom[1]-1),2) # better than Mom[1]*(Mom[1]-1)/2
    coag.fast_sample_reactants! = fast_sample_uniform_2
    coag.pi = (yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}) -> for d=1:length(yc[1]) yc[1][d] = xc[1][d] + xc[2][d] end

    frag=TransitionClass(1,2,kF)
    frag.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]
    frag.g = xc::Vector{Vector{Int64}} -> xc[1][1]
    frag.fast_sample_reactants! = fast_sample_mass_1
    frag.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                    for d=1:length(yc[1])
                        yc[1][d] = rand(0:xc[1][d])
                        yc[2][d] = xc[1][d]-yc[1][d]
                    end
              end

    intake=TransitionClass(0,1,kI)
    intake.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> 1
    intake.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}}, param::Float64)
                    for d=1:length(yc[1])
                        yc[1][d]=rand(Poisson(param))
                    end
                end
    intake.parameters = 50.0

    exit=TransitionClass(1,0,kE)
    exit.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]
    exit.g = xc::Vector{Vector{Int64}} -> 1
    exit.fast_sample_reactants! = fast_sample_uniform_1

    add_transition_class(S,coag,frag,intake,exit)

    S.moment_equations = coagulation_fragmentation_ODEs
    S.init_moments = coagulation_fragmentation_initial
    return S
end



# CASE STUDY "Transciption dynamics in a cell community"

function CellCommunity()
    S=System("CellCommunity",2,Dict(1=>[0,0],2=>[1,0],3=>[0,1]))
    kcom = 0.01
    kbG, kdG, kS, kbS, kdS = 0.01, 0.1, 1.0, 0.0, 0.05

    comm=TransitionClass(2,2,kcom)
    comm.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]*(Mom[1]-Mom[2])
    comm.g = xc::Vector{Vector{Int64}} -> xc[1][1]*(1-xc[2][1]) + xc[2][1]*(1-xc[1][1])
    comm.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                    yc[1][1], yc[2][1] = 1, 1
                    yc[1][2], yc[2][2] = xc[1][2], xc[2][2]
              end
    comm.fast_sample_reactants! = fast_sample_communication

    birthG=new_chemical_reaction_class([1,0],kbG)
    birthG.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]-Mom[2]
    birthG.g = xc::Vector{Vector{Int64}} -> 1-xc[1][1]
    #birthG.fast_sample_reactants! = fast_sample_inactive

    deathG=new_chemical_reaction_class([-1,0],kdG)
    deathG.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]
    deathG.g = xc::Vector{Vector{Int64}} -> xc[1][1]
    deathG.fast_sample_reactants! = fast_sample_mass_1

    catalysis=new_chemical_reaction_class([0,1],kS)
    catalysis.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]
    catalysis.g = xc::Vector{Vector{Int64}} -> xc[1][1]
    catalysis.fast_sample_reactants! = fast_sample_mass_1

    birthS=new_chemical_reaction_class([0,1],kbS)
    birthS.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]
    birthS.fast_sample_reactants! = fast_sample_uniform_1

    deathS=new_chemical_reaction_class([0,-1],kdS)
    deathS.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[3]
    deathS.g = xc::Vector{Vector{Int64}} -> xc[1][2]
    deathS.fast_sample_reactants! = fast_sample_mass_2

    add_transition_class(S,comm,birthG,deathG,catalysis,birthS,deathS)

    S.moment_equations = cell_community_ODEs
    S.init_moments = cell_community_initial
    return S
end



# CASE STUDY "Stem Cell Population Dynamics"

function StemCells()
    S=System("StemCells",2,Dict(1=>[0,0],2=>[1,0],3=>[1,1]))
    knf = 0.01
    kF_asym, kF_sym, kS, kEx = 0.005, 0.005, 10.0, 0.05

    nf=TransitionClass(2,2,knf)
    nf.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> div(Mom[2]*(Mom[2]-1),2) #better than Mom[2]*(Mom[2]-1)/2
    nf.g = xc::Vector{Vector{Int64}} -> xc[1][1]*xc[2][1]
    nf.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                    yc[1][1], yc[1][2] = xc[1][1], xc[1][2]
                    yc[2][1], yc[2][2] = xc[2][1], xc[2][2]
                    rand() < 0.5 ? yc[1][1]=0 : yc[2][1]=0
              end
    nf.fast_sample_reactants! = fast_sample_negative_feedback

    catalysis=new_chemical_reaction_class([0,1],kS)
    catalysis.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[2]
    catalysis.g = xc::Vector{Vector{Int64}} -> xc[1][1]
    catalysis.fast_sample_reactants! = fast_sample_mass_1

    f_asym=TransitionClass(1,2,kF_asym)
    f_asym.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[3]
    f_asym.g = xc::Vector{Vector{Int64}} -> xc[1][1]*xc[1][2]
    f_asym.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                        yc[1][1] , yc[2][1] = 1 , 0
                        yc[1][2] , yc[2][2] = 0 , 0
                end
    f_asym.fast_sample_reactants! = fast_sample_stem

    f_sym=TransitionClass(1,2,kF_sym)
    f_sym.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[3]
    f_sym.g = xc::Vector{Vector{Int64}} -> xc[1][1]*xc[1][2]
    f_sym.pi = function(yc::Vector{Vector{Int64}}, xc::Vector{Vector{Int64}})
                      yc[1][1] , yc[2][1] = 1 , 1
                      yc[1][2] , yc[2][2] = 0 , 0
                end
    f_sym.fast_sample_reactants! = fast_sample_stem

    exit=TransitionClass(1,0,kEx)
    exit.H = (n::Matrix{Int64},Mom::Vector{Int64}) -> Mom[1]-Mom[2]
    exit.g = xc::Vector{Vector{Int64}} -> 1-xc[1][1]
    exit.fast_sample_reactants!  = fast_sample_inactive

    add_transition_class(S,nf,catalysis,f_asym,f_sym,exit)

    S.moment_equations = stemcells_ODEs
    S.init_moments = stemcells_initial
    return S
end



### SUPPORTING FUNCTIONS FOR EFFICIENT SIMULATION ###


function fast_sample_uniform_1(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    r_indices[1] = rand(1:Mom[1])
end


function fast_sample_uniform_2(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    r_indices[1] = rand(1:Mom[1])
    r_indices[2] = rand(1:Mom[1])
    while r_indices[1] == r_indices[2]
        r_indices[1] = rand(1:Mom[1])
        r_indices[2] = rand(1:Mom[1])
    end
    sort!(r_indices)
end

function fast_sample_mass_1(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    rv=rand()*Mom[2]
    r_indices[1]=1
    val = 1.0*n[1,1]
    while val < rv
        r_indices[1]+=1
        val += n[1,r_indices[1]]
    end
end

function fast_sample_mass_2(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    rv=rand()*Mom[3]
    r_indices[1]=1
    val = 1.0*n[2,1]
    while val < rv
        r_indices[1]+=1
        val += n[2,r_indices[1]]
    end
end


function fast_sample_stem(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    rv=rand()*Mom[3]
    r_indices[1]=1
    val = 1.0*n[1,1]*n[2,1]
    while val < rv
        r_indices[1]+=1
        val += n[1,r_indices[1]]*n[2,r_indices[1]]
    end
end


function fast_sample_inactive(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    rv_inactive = rand(1:(Mom[1]-Mom[2]))
    r_indices[1]=1
    count_inactive = 1-n[1,1]
    while count_inactive < rv_inactive
        n[1,r_indices[1] += 1] == 0 ? count_inactive += 1 : nothing
    end
end


function fast_sample_communication(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    rv_active = rand(1:Mom[2])
    r_indices[1]=1
    count_active = n[1,1]
    while count_active < rv_active
        n[1,r_indices[1] += 1] == 1 ? count_active += 1 : nothing
    end
    rv_inactive = rand(1:(Mom[1]-Mom[2]))
    r_indices[2]=1
    count_inactive = 1-n[1,1]
    while count_inactive < rv_inactive
        n[1,r_indices[2] += 1] == 0 ? count_inactive += 1 : nothing
    end
end


function fast_sample_negative_feedback(r_indices::Vector{Int64},n::Matrix{Int64},Mom::Vector{Int64})
    r_indices[1] = rand(1:Mom[2])
    r_indices[2] = rand(1:Mom[2])
    while r_indices[1] == r_indices[2]
        r_indices[1] = rand(1:Mom[2])
        r_indices[2] = rand(1:Mom[2])
    end
    sort!(r_indices)
    cell_label = r_indices[1]
    r_indices[1]=1
    count = n[1,1]
    while count < cell_label
        n[1,r_indices[1] += 1] == 1 ? count += 1 : nothing
    end
    cell_label = r_indices[2]
    r_indices[2]=r_indices[1]
    while count < cell_label
        n[1,r_indices[2] += 1] == 1 ? count += 1 : nothing
    end
end
