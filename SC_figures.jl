
# Concept Figure

function Figure_1(;seed=5517)
seed != nothing ? Random.seed!(seed) : seed=rand(1:10^6)
S=ConceptExample()
figure_size=(2.2,2)
N0=5
n0=zeros(Int64,2,N0)
Random.seed!(seed)
for i=1:N0
	n0[1,i]=rand(Poisson(5))
end
xmax=15
ymax=12
Random.seed!(seed)
n_story,m=StochasticCompartments.SSA(S,n0,[0.;25.0;50.0],full_story=true)
for t=1:3
	figure(t,figsize=figure_size)
	HH=fit(Histogram,n_story[t][1,:],collect(0:1:xmax)); bar(collect(HH.edges[1][1:end-1]),HH.weights,width=1,alpha=0.5,color="g",label=L"x_1")
	HH=fit(Histogram,n_story[t][2,:],collect(0:1:xmax)); bar(collect(HH.edges[1][1:end-1]),HH.weights,width=1,alpha=0.5,color="orange",label=L"x_2")
    ylim([0,ymax])
    yticks([0,5,10])
end
figure(1); legend(loc="center right")
Random.seed!(seed)
n_story,t,m=StochasticCompartments.SSA(S,n0,50.0,full_story=true)
figure("compartment number",figsize=(6.5,3))
plot(t,m[1,:],drawstyle="steps-post",color="k")
xlim([0.,52.5]); ylim([0.,24.95])
figure("mass",figsize=(6.5,3))
plot(t,m[2,:],drawstyle="steps-post",color="green",label=L"x_1")
plot(t,m[3,:],drawstyle="steps-post",color="orange",label=L"x_2")
xlim([0.,52.5]); ylim([0.,97.65]); legend()
return
end


# Figures for "Nested Birth-Death Process"

function Figure_2BC(;intake_mean::Float64=10.,N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = BirthDeath()
S.c[3].parameters = intake_mean
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


function Figure_S1(;lambdas=[1.;10.;50.])
n0 = rand(0:0,1,1)
n0fix = zeros(Int64,1,100)
count=0
figure(figsize=(14,5))
for lambda in lambdas
	subplot(1,3,count+=1)  # ; title("Expected total population mass",size=14)
	S = BirthDeath(); kb=S.c[1].k ; kd=S.c[2].k
	S.c[3].parameters = lambda
	t_ode,Mom_ode=solve_moment_equations(S,n0,(0.,550.),Tpoints=1001)
	stdM=sqrt.(Mom_ode[5,:].-Mom_ode[2,:].^2)
	fill_between(t_ode, Mom_ode[2,:] .- stdM, Mom_ode[2,:] .+ stdM,alpha=0.6,color="xkcd:sky blue",label="Dynamic population")#, "*L"N(0)=1\,,\,\,\beta="*"$(lambda*kd/kb)")
	plot(t_ode, Mom_ode[2,:],color="b")
	S.c[3].k, S.c[4].k = 0. , 0.
	t_ode,Mom_ode=solve_moment_equations(S,n0fix,(0.,550.),Tpoints=1001)
	stdM=sqrt.(Mom_ode[5,:].-Mom_ode[2,:].^2)
	fill_between(t_ode, Mom_ode[2,:] .- stdM, Mom_ode[2,:] .+ stdM,alpha=0.6,color="orange",label="Fixed population")#, "*L"N(0)=100")
	plot(t_ode, Mom_ode[2,:],color="red")
	ylim([-50,1600])
	title(L"\beta="*"$(lambda*kd/kb)")
	legend(loc="lower right")
end
subplots_adjust(top=0.9,bottom=0.15,left=0.105,right=0.995,hspace=0.2,wspace=0.2)
return
end


# Figures for "Stochastic Coagulation-Fragmentation Dynamics"

function Figure_2EF(;N_SSA::Int64=10000,Coag_rates::Vector{Float64}=[0.0005;0.005;0.05],cc=["r","g","b"],seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CoagulationFragmentation();
t_ssa = [0.;5;10;15;20;40;60;80;100]
n0 = rand(10:10,1,100)
kF=S.c[2].k
for j=3:-1:1
    S.c[1].k=Coag_rates[j]
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
figure(2); yticks([0,1000,2000,3000,4000,5000]); ylim([-100,6150]); legend();
figure(3); ax=gca(); ax.set_yscale("log")
return
end



# Figures for "Transciption dynamics in a cell community"

function Figure_3B(;
    Comm_rates=[0.;0.001;0.002;0.005;0.01;0.05],
    farbe=["grey","orange","red","green","blue","purple"],
    N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CellCommunity()
N_0 = 100
t_ssa = [0.;5;10;15;20;30;40;60;80;100;125;150]
n0 = zeros(Int64,2,N_0); n0[:,1] .= 1
figure(1); title("Expected number of active cells",size=14); #xlabel("time",size=12)
figure(2); title("Expected total protein amount ",size=14); #xlabel("time",size=12)
kdG = S.c[3].k
for index=length(Comm_rates):-1:1
    S.c[1].k = Comm_rates[index]
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
figure(1); xticks([0,25,50,75,100,125,150]); ylim([-5,105]); subplots_adjust(top=0.9,bottom=0.11,left=0.11,right=0.9,hspace=0.2,wspace=0.2)
figure(2); xticks([0,25,50,75,100,125,150]); yticks([0,500,1000,1500,2000]); ylim([-100,2100]); subplots_adjust(top=0.9,bottom=0.11,left=0.11,right=0.9,hspace=0.2,wspace=0.2)
return
end


function Figure_3C(;
    time_interval = [0.;200.],
    Comm_rates=10.0 .^ collect(-5:0.01:1),
    seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = CellCommunity()
n0 = zeros(Int64,2,100); n0[1,1] = 1
MGs=zeros(length(Comm_rates))
stdMGs=zeros(length(Comm_rates))
MPs=zeros(length(Comm_rates))
stdMPs=zeros(length(Comm_rates))
for index=1:length(Comm_rates)
    S.c[1].k = Comm_rates[index]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(time_interval[1],time_interval[end]),Tpoints=1001)
	MGs[index]=Mom_ode[2,end]
    stdMGs[index]=sqrt.(Mom_ode[4,:].-Mom_ode[2,:].^2)[end]
	MPs[index]=Mom_ode[3,end]
    stdMPs[index]=sqrt.(Mom_ode[5,:].-Mom_ode[3,:].^2)[end]
end
Comm_rescaled=size(n0,2)*Comm_rates/S.c[3].k
# active cells
figure(); xscale("log"); #xlabel(L"\frac{k_{com} N_0}{k_d^G}",size=14); title("Expected steady-state number of active cells",size=14)
fill_between(Comm_rescaled, MGs .- stdMGs, MGs .+ stdMGs,alpha=0.5,color="xkcd:sky blue")
plot(Comm_rescaled, MGs, color="xkcd:royal blue")
ylim([-5,105])
yticks([0,20,40,60,80,100])
ax=gca()
ax2 = ax.twinx() # Create another axis on top of the current axis
setp(ax.get_yticklabels(),color="xkcd:royal blue") # Y Axis font formatting
plot(Comm_rescaled, stdMGs.^2 ./ MGs, "-.", lw=1.0, color="k")
yticks([0,0.25,0.5,0.75,1.,1.25,1.5]); #legend(loc="center right",fontsize=16)
subplots_adjust(top=0.9,bottom=0.11,left=0.11,right=0.9,hspace=0.2,wspace=0.2)
# protein mass
figure(); xscale("log"); #xlabel(L"\frac{k_{com} N_0}{k_d^G}",size=14); title("Expected steady-state protein mass",size=14)
ax=gca()
plot(Comm_rescaled, stdMPs.^2 ./ MPs, "-", lw=1.0, color="k")#, label=L"\frac{\sigma^2}{\mu}")
xlim([0.01,2500])
yticks([0,1,5,10,15]); #legend(loc="center right",fontsize=16)
subplots_adjust(top=0.9,bottom=0.11,left=0.11,right=0.9,hspace=0.2,wspace=0.2)
return
end


function Figure_3C_inset(; tmax::Float64=150.,
                                Nsamples::Int64=1000,
								edges=collect(0.:1.:100),
								Comm_rates=[0.;0.001;0.002;0.005;0.01;0.05],
							    farbe=["grey","orange","red","green","blue","purple"],
                                max_allocation::Int64=10^6,
								color="green",
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
tt=[0.;tmax]
n0 = zeros(Int64,2,100); n0[:,1] .= 1
S = CellCommunity()
kb1,kd1,kp,kd2 = [S.c[i].k for i in[2;3;4;6] ]
dh = edges[2]-edges[1]
h_centers = edges[1:end-1]
figure("Steady-state protein distribution")
for i=1:length(Comm_rates)
	S.c[1].k = Comm_rates[i]
	kcomm_eff = Comm_rates[i]*100/S.c[3].k
	samples=zeros(Int64,max_allocation)
	writing=0
	for i=1:Nsamples
	    println(i, " / ",Nsamples)
	    n,MM=SSA(S,n0,tt)
	    for i=1:size(n,2)
			writing +=1
			writing > length(samples) ? samples=[samples;zeros(Int64,max_allocation)] : nothing
			samples[writing]=n[2,i]
	    end
	end
	samples=samples[1:writing]
	HH=fit(Histogram,samples,edges)
	histo=HH.weights/dh/sum(HH.weights)
	plot(h_centers,HH.weights/dh/sum(HH.weights),alpha=0.5,color=farbe[i],label="$kcomm_eff")
end
plot(h_centers,pdf_transcription_telegraph.(h_centers,kb1,kd1,kp,kd2),":",color="k",label="analytical, "*L"k_{com}=0")
plot(h_centers,pdf.(Poisson(kp/kd2),h_centers),"--",color="k",label="analytical, "*L"k_{com} = \infty")
xlabel(L"x_S",size=12); xlim([-0.5,35.5]); ylim([-0.01,0.31]); ylabel(L"P(x_S)",size=12,rotation=90)
#legend()
return
end

function pdf_transcription_telegraph(x,kb1,kd1,kp,kd2)
b1d2, d1d2, b2d2 = kb1/kd2, kd1/kd2, kp/kd2
log_gammas = loggamma(b1d2+x)-loggamma(x+1)-loggamma(b1d2+d1d2+x)+loggamma(b1d2+d1d2)-loggamma(b1d2)
return exp(log_gammas)*(b2d2)^x*mFn([b1d2+x],[b1d2+d1d2+x],-b2d2)
end


# Figures for "Stem Cell Population Dynamics"

function Figure_3E(; tmax::Float64=20.,
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


function Figure_3F(;N_SSA::Int64=1000,seed::Union{Nothing,Int64}=nothing,show_leg::Bool=true)
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


function Figure_3H(;seed::Union{Nothing,Int64}=75398,show_leg::Bool=true) # 5734 # 2517 # 8512
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
n0 = ones(Int64,2,1); n0[1]=1
figure("f_ST(t)",figsize=(5,4))
t_ssa1=collect(0:0.1:200.)
n1,m1=StochasticCompartments.SSA(S,n0,t_ssa1)
t_ode1,Mom_ode1=solve_moment_equations(S,n0,(t_ssa1[1],t_ssa1[end]),Tpoints=1001)
t_ssa2 = t_ssa1 .+ 200.
S.c[1].k /= 5
n2,m2=StochasticCompartments.SSA(S,n1,t_ssa2)
t_ode2,Mom_ode2=solve_moment_equations(S,Mom_ode1[:,end],(t_ssa2[1],t_ssa2[end]),Tpoints=1001)
tt = [t_ssa1;t_ssa2]
mm = [m1 m2]
plot(tt,mm[2,:] ./ mm[1,:], lw=1., color="orange",label="SSA sim.")
t_ode = [t_ode1;t_ode2]
Mom_ode = [Mom_ode1 Mom_ode2]
f_ST = Mom_ode[2,:] ./ Mom_ode[1,:]
plot(t_ode, f_ST, lw=1.5, color="red",label="ODEs")
plot([0;350],f_ST[end]*ones(2),lw=1.,"--",color="k")
xlim([135,365])
xticks([150,200,250,300,350])
ylim([0.05,0.45])
yticks([0.1,0.2,0.3,0.4])
legend()
return
end


function Figure_3G(;ks0::Float64=10.,knf0::Float64=0.01,theta0::Float64=1.,
                                changing::Vector{Float64}=[0.1;0.2;0.5;1;2;5;10],
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.c[1].k = knf0
S.c[2].k = ks0
KF = S.c[3].k + S.c[4].k
S.c[3].k = KF/(1+theta0)
S.c[4].k = theta0*S.c[3].k
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
S.c[1].k = knf0
S.c[2].k = ks0
KF = S.c[3].k + S.c[4].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(theta)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
for i=1:LL
	S.c[3].k = KF/(1+theta[i])
	S.c[4].k = theta[i]*S.c[3].k
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
return Ntot, Stot, Nstem, Sstem
end



function StemCells_figure_ks(; theta0::Float64=1., knf0::Float64=0.01,
                            ks::Vector{Float64}=[1.;2.;3;5.;7.;10;15;20;30;50;70;100],
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.c[1].k = knf0
KF = S.c[3].k + S.c[4].k
S.c[3].k = KF/(1+theta0)
S.c[4].k = theta0*S.c[3].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(ks)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
for i=1:LL
    S.c[2].k = ks[i]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
return Ntot, Stot, Nstem, Sstem
end



function StemCells_figure_knf(;theta0::Float64=1.,ks0::Float64=10.,
                              knf::Vector{Float64}= 10.0 .^ collect(-3:0.25:-1),
                                seed::Union{Nothing,Int64}=nothing)
seed!=nothing ? Random.seed!(seed) : nothing
S = StemCells()
S.c[2].k = ks0
KF = S.c[3].k + S.c[4].k
S.c[3].k = KF/(1+theta0)
S.c[4].k = theta0*S.c[3].k
t_ssa = collect(0.:50:1000)
n0 = ones(Int64,2,1)
LL=length(knf)
Nstem=zeros(LL); Sstem=zeros(LL); Ntot=zeros(LL); Stot=zeros(LL)
for i=1:LL
    S.c[1].k = knf[i]
    t_ode,Mom_ode=solve_moment_equations(S,n0,(t_ssa[1],t_ssa[end]),Tpoints=1001)
    Ntot[i]=Mom_ode[1,end]
    Stot[i]=sqrt.(Mom_ode[5,end].-Mom_ode[1,end].^2)
    Nstem[i]=Mom_ode[2,end]
    Sstem[i]=sqrt.(Mom_ode[6,end].-Mom_ode[2,end].^2)
end
return Ntot, Stot, Nstem, Sstem
end



function Figure_S4(; kb=10.,kF=0.01,t=collect(0:0.1:15),Nsamples=1000000)
taus=zeros(Nsamples)
for i=1:Nsamples
	time=0.
	x=0
	while time<t[end]
		h_tot = kb+kF*x
		wt = -log(1-rand())/h_tot
		if rand() < kb/h_tot
			time+= wt
			x+=1
		else
			taus[i]=time+wt
			break
		end
	end
end
#plot(t,kF*kb*t .* exp.(-kF*kb*t.^2 / 2),color="k","--",label="approximated analytical result")
HH=fit(Histogram,taus,t)
dh=HH.edges[1][2]-HH.edges[1][1]
centers=collect(HH.edges[1][1:end-1])
histo=HH.weights/dh/sum(HH.weights)
bar(collect(HH.edges[1][1:end-1]),HH.weights/dh/sum(HH.weights),width=dh,alpha=0.5,color="b")#,label="sampled histogram")
xlabel(L"\tau",rotation=0,size=16,labelpad=5)
ylabel(L"p(\tau)",rotation=0,size=16,labelpad=25)
title("Inter-division-time distribution",size=16)
return
end
