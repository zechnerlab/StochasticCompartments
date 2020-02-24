
# Concept Figure

function figure_1(;rg_seed::Int64=5517)
S=ConceptFigure()
figure_size=(2.2,2)
N0=100
n0=zeros(Int64,2,N0)
Random.seed!(rg_seed)
for i=1:N0
	n0[1,i]=rand(Poisson(5))
	n0[2,i]=rand(Poisson(2))
end
xmax=30
ymax=30
Random.seed!(rg_seed)
n_story,m=StochasticCompartments.SSA(S,n0,[0.;1.0;2.0],full_story=true)
for t=1:3
	figure(t,figsize=figure_size)
	HH=fit(Histogram,n_story[t][1,:],collect(0:1:xmax)); bar(collect(HH.edges[1][1:end-1]),HH.weights,width=1,alpha=0.5,color="g",label=L"x_1")
	HH=fit(Histogram,n_story[t][2,:],collect(0:1:xmax)); bar(collect(HH.edges[1][1:end-1]),HH.weights,width=1,alpha=0.5,color="orange",label=L"x_2")
    ylim([0,ymax])
    yticks([0,5,10,15,20,25])
end
figure(1); legend(loc="center right")
return
end


# Figures for "Nested Birth-Death Process"

function figure_birthdeath(;intake_mean::Float64=10.,N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = BirthDeath()
S.transition_classes[3].parameters = intake_mean
t_ssa = collect(0.:25:500)
n0 = rand(1:1,1,1)
Mom_ssa,Var_ssa=SSA(S,n0,t_ssa,N_SSA)
t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
stdN=sqrt.(Mom_ode[4,:].-Mom_ode[1,:].^2)
stdM=sqrt.(Mom_ode[5,:].-Mom_ode[2,:].^2)
for i=1:2
    figure(i)
    errorbar(t_ssa,Mom_ssa[i,:], yerr=[sqrt.(Var_ssa[i,:]) sqrt.(Var_ssa[i,:])]',fmt="o",ms=5.,capsize=5.,color="blue",label="SSA",alpha=0.5)
end
figure(1); title("Average compartment number",size=14)
fill_between(t_ode, Mom_ode[1,:] .- stdN, Mom_ode[1,:] .+ stdN,alpha=0.6,color="xkcd:sky blue",label="ODEs")
plot(t_ode, Mom_ode[1,:],color="b")
figure(2); title("Expected total population mass",size=14)
fill_between(t_ode, Mom_ode[2,:] .- stdM, Mom_ode[2,:] .+ stdM,alpha=0.6,color="xkcd:sky blue",label="ODEs")
plot(t_ode, Mom_ode[2,:],color="b")
for i=1:2 figure(i); legend(); xlabel("time",size=12) end
return (t_ssa,Mom_ssa,Var_ssa), (t_ode,Mom_ode)
end


function figure_SI1(;beta::Vector{Float64}=[1.;2;5;7;10;15;20;30;40;50;70;100]/10,alpha=[0.01; 0.05; 0.1; 0.5; 1.;])
figure("VMR")
S = BirthDeath()
kb = S.transition_classes[1].k
kd = S.transition_classes[2].k
for i=1:length(alpha)
    kE = alpha[i]*kd
    factor = kb/kd*alpha[i]/(1+alpha[i])/(2+alpha[i])
    plot(beta, 1 .+ factor*(beta .- 1).^2 ./ (1 .+ beta*alpha[i]).^2 ,label=L"\alpha" * " = $(alpha[i])")
end
xscale("log"); xticks([0.1,1,10],["0.1","1","10"],size=12); xlim([0.07;12])
yscale("log"); yticks([1,10],["1","10"]); ylim([0.7,12])
xlabel("Rescaled intake mean, " * L"\beta",size=12)
ylabel(L"\frac{\mathrm{Var}(X_\infty)}{\mathbb{E}[X_\infty]}",rotation=0,size=20,labelpad=15)
subplots_adjust(top=0.88,bottom=0.11,left=0.14,right=0.93,hspace=0.2,wspace=0.2)
title("Variance-to-mean ratio of compartment content at steady-state",size=12)
legend()
return
end



# Figures for "Stochastic Coagulation-Fragmentation Dynamics"

function figure_coagulationfragmentation(;N_SSA::Int64=10000,Coag_rates::Vector{Float64}=[0.0005;0.005;0.05],cc=["r","g","b"],seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CoagulationFragmentation();
t_ssa = [0.;5;10;15;20;40;60;80;100]
n0 = rand(10:10,1,100)
kF=S.transition_classes[2].k
for j=3:-1:1
    S.transition_classes[1].k=Coag_rates[j]
    cond_name=L"k_C/k_F=" * "$(Coag_rates[j]/kF)"
    Mom_ssa,Var_ssa=SSA(S,n0,t_ssa,N_SSA)
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    stdN=sqrt.(Mom_ode[4,:].-Mom_ode[1,:].^2)
    stdM=sqrt.(Mom_ode[5,:].-Mom_ode[2,:].^2)
    for i=1:3
        figure(i)
        errorbar(t_ssa,Mom_ssa[i,:], yerr=[sqrt.(Var_ssa[i,:]) sqrt.(Var_ssa[i,:])]',fmt="o",ms=5.,capsize=5.,color=cc[j],label=cond_name * " (SSA)",alpha=0.9)
    end
    figure(1); title("Average compartment number",size=12)
    fill_between(t_ode, Mom_ode[1,:] .- stdN, Mom_ode[1,:] .+ stdN,alpha=0.4,color=cc[j],label=cond_name * " (ODEs)")
    plot(t_ode, Mom_ode[1,:],color=cc[j])
    figure(2); title("Expected total population mass",size=12)
    fill_between(t_ode, Mom_ode[2,:] .- stdM, Mom_ode[2,:] .+ stdM,alpha=0.4,color=cc[j],label=cond_name * " (ODEs)")
    plot(t_ode, Mom_ode[2,:],color=cc[j])
    figure(3); title("Expected total squared content",size=12)
    plot(t_ode, Mom_ode[3,:],color=cc[j],label=cond_name * " (ODEs)")
end
for i=1:3 figure(i); xlabel("time",size=12) end
figure(1); ylim([-5,265])
figure(2); legend(); yticks([1000,2000,3000,4000,5000])
figure(3); ax=gca(); ax.set_yscale("log")
return
end



# Figures for "Transciption dynamics in a cell community"

function figure_cellcommunication(;
    Comm_rates=[0.;0.001;0.002;0.005;0.01;0.05],
    farbe=["grey","orange","red","green","blue","purple"],
    N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CellCommunity();
N_0 = 100
t_ssa = [0.;5;10;15;20;30;40;60;80;100;125;150]
n0 = zeros(Int64,2,N_0); n0[1,1] = 1
figure(1); title("Expected number of active cells",size=14)
figure(2); title("Expected total protein amount ",size=14)
kdG = S.transition_classes[3].k
for index=length(Comm_rates):-1:1
    S.transition_classes[1].k = Comm_rates[index]
    Mom_ssa,Var_ssa=SSA(S,n0,t_ssa,N_SSA)
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    stdMG=sqrt.(Mom_ode[4,:].-Mom_ode[2,:].^2)
    stdMP=sqrt.(Mom_ode[5,:].-Mom_ode[3,:].^2)
    figure(1)
    errorbar(t_ssa,Mom_ssa[2,:], yerr=[sqrt.(Var_ssa[2,:]) sqrt.(Var_ssa[2,:])]',fmt="o",ms=5.,capsize=5.,color=farbe[index])
    fill_between(t_ode, Mom_ode[2,:] .- stdMG, Mom_ode[2,:] .+ stdMG,alpha=0.5,color=farbe[index])
    plot(t_ode, Mom_ode[2,:],color=farbe[index],label=L"k_{com} N_0/k_d^G =" * "$(round(Int,Comm_rates[index]*N_0/kdG))")
    figure(2)
    errorbar(t_ssa,Mom_ssa[3,:], yerr=[sqrt.(Var_ssa[3,:]) sqrt.(Var_ssa[3,:])]',fmt="o",ms=5.,capsize=5.,color=farbe[index])
    fill_between(t_ode, Mom_ode[3,:] .- stdMP, Mom_ode[3,:] .+ stdMP,alpha=0.5,color=farbe[index])
    plot(t_ode, Mom_ode[3,:],color=farbe[index],label=L"k_{com} N_0/k_d^G =" * "$(round(Int,Comm_rates[index]*N_0/kdG))")
end
for i=1:2 figure(i); legend(loc="lower right"); xlabel("time",size=12) end
return
end


function figure_cellactivation(;
    time_interval = [0.;200.],
    #Comm_rates=[0.00001;0.00002;0.00005;0.0001;0.0002;0.0005;0.001;0.002;0.005;0.01;0.02;0.05;0.1;0.2;0.5;1.],
    Comm_rates=10.0 .^ collect(-5:0.1:0),
    seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CellCommunity();
n0 = zeros(Int64,2,100); n0[1,1] = 1
MGs=zeros(length(Comm_rates))
stdMGs=zeros(length(Comm_rates))
for index=1:length(Comm_rates)
    S.transition_classes[1].k = Comm_rates[index]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(time_interval[1],time_interval[end]),Tpoints=1001)
    MGs[index]=Mom_ode[2,end]
    stdMGs[index]=sqrt.(Mom_ode[4,:].-Mom_ode[2,:].^2)[end]
end
Comm_rescaled=size(n0,2)*Comm_rates/S.transition_classes[3].k
figure(); xscale("log"); xlabel(L"\frac{k_{com} N_0}{k_d^G}",size=14); title("Expected steady-state number of active cells",size=14)
fill_between(Comm_rescaled, MGs .- stdMGs, MGs .+ stdMGs,alpha=0.5,color="xkcd:sky blue")
plot(Comm_rescaled, MGs, color="xkcd:royal blue")
yticks([0,20,40,60,80,100])
ax=gca()
ax2 = ax.twinx() # Create another axis on top of the current axis
setp(ax.get_yticklabels(),color="xkcd:royal blue") # Y Axis font formatting
plot(Comm_rescaled, stdMGs.^2 ./ MGs, "-.", lw=1.0, color="k")#, label=L"\frac{\sigma^2}{\mu}")
#ylabel(L"\frac{\sigma^2}{\mu}",size=16, color="xkcd:royal blue",rotation=0,labelpad=10)
yticks([0,0.25,0.5,0.75,1.,1.25,1.5]); #legend(loc="center right",fontsize=16)
subplots_adjust(top=0.9,bottom=0.155,left=0.11,right=0.9,hspace=0.2,wspace=0.2)
return
end



# Figures for "Stem Cell Population Dynamics"

function figure_stemcelldynamics(;N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing,show_leg::Bool=true)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells();
t_ssa = collect(0.:12.5:200)
n0 = ones(Int64,2,1); n0[1]=1
# n0 = ones(Int64,2,100)
Mom_ssa,Var_ssa=SSA(S,n0,t_ssa,N_SSA)
t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
stdN=sqrt.(Mom_ode[5,:].-Mom_ode[1,:].^2)
stdMG=sqrt.(Mom_ode[6,:].-Mom_ode[2,:].^2)
figure("Stem-vs-total")
errorbar(t_ssa,Mom_ssa[1,:], yerr=[sqrt.(Var_ssa[1,:]) sqrt.(Var_ssa[1,:])]',fmt="o",ms=5.,capsize=5.,color="b",label="All cells - SSA",alpha=0.5)
plot(t_ode, Mom_ode[1,:], color="blue")
fill_between(t_ode, Mom_ode[1,:] .- stdN, Mom_ode[1,:] .+ stdN,alpha=0.6,color="xkcd:sky blue",label="All cells - ODEs")
errorbar(t_ssa,Mom_ssa[2,:], yerr=[sqrt.(Var_ssa[2,:]) sqrt.(Var_ssa[2,:])]',fmt="o",ms=5.,capsize=5.,color="red",label="Stem cells - SSA",alpha=0.5)
plot(t_ode, Mom_ode[2,:], color="red")
fill_between(t_ode, Mom_ode[2,:] .- stdMG, Mom_ode[2,:] .+ stdMG,alpha=0.6,color="xkcd:orange",label="Stem cells - ODEs")
show_leg ? legend() : nothing
xlabel("time",size=12)
return (t_ssa,Mom_ssa,Var_ssa), (t_ode,Mom_ode)
end


function figure_stemcellparameters(;ks0::Float64=10.,knf0::Float64=0.01,theta0::Float64=1.,
                                changing::Vector{Float64}=[0.1;0.2;0.5;1;2;5;10],
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.transition_classes[1].k = knf0
S.transition_classes[2].k = ks0
KF = S.transition_classes[3].k + S.transition_classes[4].k
S.transition_classes[3].k = KF/(1+theta0)
S.transition_classes[4].k = theta0*S.transition_classes[3].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
nt_theta,st_theta,ns_theta,ss_theta=StemCells_figure_theta(theta=changing*theta0,knf0=knf0,ks0=ks0)
fraction_theta = ns_theta ./ nt_theta
nt_ks,st_ks,ns_ks,ss_ks=StemCells_figure_ks(ks=changing*ks0,knf0=knf0,theta0=theta0)
fraction_ks = ns_ks ./ nt_ks
nt_nf,st_nf,ns_nf,ss_nf=StemCells_figure_knf(knf=changing*knf0,ks0=ks0,theta0=theta0)
fraction_nf = ns_nf ./ nt_nf
#FRACTION#
figure("parameter_sweep_fraction",figsize=(5,4)); title("Steady-state fraction of stem cells")
ax=gca(); ax.set_xscale("log"); ylim([0,0.5])
plot(changing, fraction_theta, "D-", label="varying "*L"\theta^+_{-}",color="blue")
plot(changing, fraction_nf, "D-", label="varying "*L"k_{com}",color="green")
plot(changing, fraction_ks, "D-", label="varying "*L"k_S",color="red")
legend(loc="upper right"); xlabel("Multypling factor", size=12)
#ABSOLUTE#
figure("parameter_sweep_absolute",figsize=(5,4)); title("Steady-state number of stem cells")
ax=gca(); ax.set_xscale("log"); ylim([1,300]); ax.set_yscale("log")
plot(changing, ns_theta, color="blue", "D-",label="varying "*L"\theta^+_{-}")
fill_between(changing, ns_theta .- ss_theta, ns_theta .+ ss_theta,alpha=0.3,color="blue")
plot(changing, ns_nf, color="green","D-",label="varying "*L"k_{com}")
fill_between(changing, ns_nf .- ss_nf, ns_nf .+ ss_nf,alpha=0.3,color="green")
plot(changing, ns_ks, color="red","D-",label="varying "*L"k_S")
fill_between(changing, ns_ks .- ss_ks, ns_ks .+ ss_ks,alpha=0.3,color="red")
legend(loc="upper right"); xlabel("Multypling factor", size=12)
return changing, fraction_theta, fraction_ks, fraction_nf
end



function StemCells_figure_theta(; ks0::Float64=10.,knf0::Float64=0.01,
                                theta::Vector{Float64}=[0.01;0.05;0.1;0.15;0.2;0.25;0.3;0.4;0.5;0.6;0.7;0.75;0.8;0.85;0.9;0.95;0.99],
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.transition_classes[1].k = knf0
S.transition_classes[2].k = ks0
KF = S.transition_classes[3].k + S.transition_classes[4].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(theta)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
#figure(1)
for i=1:LL
	S.transition_classes[3].k = KF/(1+theta[i])
	S.transition_classes[4].k = theta[i]*S.transition_classes[3].k
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    #plot(t_ode,Mom_ode[1,:],color="b")
    #plot(t_ode,Mom_ode[2,:],color="orange")
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
#figure(2)
#plot(theta, Ntot, color="blue", "D-")
#fill_between(theta,Ntot .- Stot, Ntot .+ Stot,alpha=0.6,color="xkcd:sky blue",label="All cells - ODEs")
#plot(theta, Nstem, color="red", "D-")
#fill_between(theta, Nstem .- Sstem, Nstem .+ Sstem,alpha=0.6,color="xkcd:orange",label="Stem cells - ODEs")
#legend(); xlabel(L"\theta",size=12)
#figure(3); title("Stem cell fraction")
#plot(theta, Nstem ./ Ntot, "D-")
#xlabel(L"\theta",size=12); ylim([0,0.5])
return Ntot, Stot, Nstem, Sstem
end



function StemCells_figure_ks(; theta0::Float64=1., knf0::Float64=0.01,
                            ks::Vector{Float64}=[1.;2.;3;5.;7.;10;15;20;30;50;70;100],
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.transition_classes[1].k = knf0
KF = S.transition_classes[3].k + S.transition_classes[4].k
S.transition_classes[3].k = KF/(1+theta0)
S.transition_classes[4].k = theta0*S.transition_classes[3].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(ks)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
#figure(1)
for i=1:LL
    S.transition_classes[2].k = ks[i]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    #plot(t_ode,Mom_ode[1,:],color="b")
    #plot(t_ode,Mom_ode[2,:],color="b")
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
#figure(2)
#plot(ks, Ntot, color="blue","D-")
#fill_between(ks, Ntot .- Stot, Ntot .+ Stot,alpha=0.6,color="xkcd:sky blue",label="All cells - ODEs")
#plot(ks, Nstem, color="red","D-")
#fill_between(ks, Nstem .- Sstem, Nstem .+ Sstem,alpha=0.6,color="xkcd:orange",label="Stem cells - ODEs")
#legend(); xlabel(L"k_S",size=12)
#figure(3); title("Stem cell fraction")
#plot(ks, Nstem ./ Ntot, "D-")
#xlabel(L"k_S",size=12); ylim([0,0.5])
return Ntot, Stot, Nstem, Sstem
end



function StemCells_figure_knf(;theta0::Float64=1.,ks0::Float64=10.,
                              knf::Vector{Float64}= 10.0 .^ collect(-3:0.25:-1),
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.transition_classes[2].k = ks0
KF = S.transition_classes[3].k + S.transition_classes[4].k
S.transition_classes[3].k = KF/(1+theta0)
S.transition_classes[4].k = theta0*S.transition_classes[3].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(knf)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
#figure(1)
for i=1:LL
    S.transition_classes[1].k = knf[i]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    #plot(t_ode,Mom_ode[1,:],color="b")
    #plot(t_ode,Mom_ode[2,:],color="b")
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
#figure(2)
#plot(knf, Ntot, color="blue", "D-")
#fill_between(knf,Ntot .- Stot, Ntot .+ Stot,alpha=0.6,color="xkcd:sky blue",label="All cells - ODEs")
#plot(knf, Nstem, color="red","D-")
#fill_between(knf, Nstem .- Sstem, Nstem .+ Sstem,alpha=0.6,color="xkcd:orange",label="Stem cells - ODEs")
#legend(); ax=gca(); ax.set_xscale("log"); ax.set_yscale("log")
#xlabel(L"k_{nf}",size=12)
#figure(3); title("Stem cell fraction")
#plot(knf, Nstem ./ Ntot, "D-")
#ax=gca(); ax.set_xscale("log"); ylim([0,0.5])
#xlabel(L"k_{nf}",size=12)
return Ntot, Stot, Nstem, Sstem
end



function figure_stemcellstart(; tmax::Float64=20.,
                           seed::Union{Nothing,Int64}=nothing) # 98455 81171 60795 11091
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
n0=ones(Int64,2,1)
n_story,t,mom= StochasticCompartments.SSA(S,n0,20.,full_story=true)
figure(figsize=(7,4))
subplot(2,1,1)
for cell_id=1:maximum(mom[1,:])
    enough_cells_times=findall(mom[1,:] .>= cell_id)
	if length(enough_cells_times) > 1
		xs = zeros(Int64,length(t))
		for i in enough_cells_times
			xs[i] = n_story[i][1,cell_id] .* n_story[i][2,cell_id]
		end
		plot(t,xs,drawstyle="steps-post",color="g")
    end
end
ax1 = gca()
ylabel(L"x_S"*" in stem cells",size=13)
ylim([0,maximum(maximum.(n_story))+10])
subplot(2,1,2,sharex=ax1)
plot(t,mom[1,:],drawstyle="steps-post",label=L"N")
plot(t,mom[2,:],drawstyle="steps-post",label=L"M^{1,0}")
ylim([0,maximum(mom[1,:])+1]) # ; yticks([0,1,2,3,4,5])
legend()
xlabel("time",size=13)
subplots_adjust(top=0.9,bottom=0.12,left=0.11,right=0.92,hspace=0.0,wspace=0.2)
return
end



function figure_stemcellperturbation(;seed::Union{Nothing,Int64}=75398,show_leg::Bool=true) # 5734 # 2517 # 8512
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.moment_equations = stemcells_ODEs_perturbation
t_ssa = collect(0.:12.5:400)
n0 = ones(Int64,2,1); n0[1]=1
figure("f_ST(t)",figsize=(5,4))
t_ssa1=collect(0:0.1:200.)
n1,m1=StochasticCompartments.SSA(S,n0,t_ssa1)
t_ssa2 = t_ssa1 .+ 200.
S.transition_classes[1].k /= 5
n2,m2=StochasticCompartments.SSA(S,n1,t_ssa2)
tt = [t_ssa1;t_ssa2]
mm = [m1 m2]
plot(tt,mm[2,:] ./ mm[1,:], lw=1., color="orange",label="SSA sim.")
t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
f_ST = Mom_ode[2,:] ./ Mom_ode[1,:]
plot(t_ode, f_ST, lw=1.5, color="red",label="ODEs")
plot([0;350],f_ST[500]*ones(2),lw=1.,"--",color="k")
xlim([135,365])
xticks([150,200,250,300,350])
ylim([0.05,0.45])
yticks([0.1,0.2,0.3,0.4])
legend()
return
end


