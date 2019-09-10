module vscantom
###############################################################
# Copyright 2019 Rowan Iskandar
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#################################################################
#vscantom module include the following modeling methods to represent

##################################################################
#use packages
using DifferentialEquations, DiffEqBiological
using LinearAlgebra, Statistics, Compat
using Distributions #, Printf,Random
#using LightGraphs, MetaGraphs, GraphPlot, Colors, Compose
# using Fontconfig, Cairo
using PoissonRandom
using Dates
#SetParams initialize module with required global parameters
#t0:pop_initial time, tfin:final time, dt:time step
#dt_scale:scale dt for SDE
#c: vector of rate constants, dim(c): number of transitions
#pop_init: vector of pop_initial state-config, dim(pop_init): number of states
#nMC: number of monte-carlo simulations
#d: matrix of stoichiometries, dim(d): dim(c)x s
#N: N(t) current state-config, or N(t)
#Q: generator (Q) matrix which contains elements of c, must be square
#################################################################
function SetParamsBasic(state_names_in,t0_in,tfin_in,dt_in,age_start_in, pop_init_in,popsize_in,rDie_in)
     global state_names=state_names_in
     global t0 = t0_in
     global tfin=tfin_in
     global dt=dt_in #virus timescale
     global age_start = age_start_in #starting age
     global pop_init=pop_init_in
     global popsize=popsize_in
     global Q=Q_in
     global num_s = length(state_names)
     global tpoints = collect(t0:dt:tfin)
     global num_t = length(tpoints)
     global index_pull=collect(1:dt*dt_scale:tfin*dt_scale+dt_scale-1)
     global rDie=rDie_in
end

#parameters virus
function SetParamsVirus(beta_in,sigma_in,pIR_in,pID_in, pIncVirusDis_in, pDeathVirusDis_in)
     global beta=beta_in #infectivity
     global sigma=sigma_in
     global pIR=pIR_in
     global pID=pID_in
     global pIncVirusDis=pIncVirusDis_in
     global pDeathVirusDis=pDeathVirusDis_in
end

#parameters cancer
function SetParamsCancer(rIncCancer_in, rIncCancerRR_in, rSoujournPre_in, rSoujourn_in, pStage_in, pStageScreen_in,
    rDieCancerStage1_in, rDieCancerStage2_in, rDieCancerStage3_in, rDieCancerStage4_in)
     global rIncCancer = rIncCancer_in
     global rIncCancerVirus = rIncCancer*rIncCancerRR_in
     global rSojournPre = rSojournPre_in
     global rSojourn = rSoujourn_in
     global pStage = pStage_in
     global pStageScreen = pStageScreen_in
     global rDieCancerStage1 = rDieCancerStage1_in
     global rDieCancerStage2 = rDieCancerStage2_in
     global rDieCancerStage3 = rDieCancerStage3_in
     global rDieCancerStage4 = rDieCancerStage4_in
end

function SetParamsScreen(pScreenSens_in, pScreenSpec_in)
     global pScreenSens = pScreenSens_in
     global pScreenSpec = pScreenSpec_in
     global intScreen = intScreen_in
end
#################################################################
function simulate(nsims_in)
#outputs dim: nsims x times
    nsims=nsims_in
    ntimes=tfin-t0+1
    n=popsize
    for l=1:nsims #monte carlo samples loop
        #time stamps for event for each individual
        #colums: 1: tDeath, 2: tIncVirus, 3: tIncCancer, 4: tDiagCancer, 5: tIncVirusDis, 6: tDeathCancer, 7:tDeathVirusDis
        tstamps = zeros(n^2,7)
        #current (clock) time
        tnow = t0
        age_now = age_start
        #setting up neighboring states matrix
        I = reshape(1:n^2,n,n)
        Iext = [NaN for i=1:n+2 for j=1:n+2]
        Iext[2:n+1,2:n+1]=I
        Neigh=[NaN for i=1:n^2 for j=1:4]
        #west neighbors
        aux=Iext[2:n+1,1:n]
        Neigh[:,1]=aux[:]
        #east neibhgors
        aux=Iext[2:n+1,3:n+2]
        Neigh[:,2]=aux[:]
        #south neibhgors
        aux=Iext[3:n+2,2:n+1]
        Neigh[:,3]=aux[:]
        #north neibhgors
        aux=Iext[1:n,2:n+1]
        Neigh[:,4]=aux[:]
        #initial distribution of infected
        ninfect=10
        #state vectors for virus part
        X=ones[n^2,1]
        X[1:2]=2; X[end]=2
        X[round(n^2/2)]=2
        X[round(n^2/3)]=2
        X[round(n^2/4)]=2
        X[round(n^2/5)]=2
        X[round(n^2/6)]=2
        X[round(n^2/7)]=2
        X[round(n^2/8)]=2
        Xnext=NaN(n^2,1)
        Xnext=[NaN for i=1:n^2 for j=1:1]
        #state vectors for cancer part: individuals*ntimes
        Y=ones[n^2,1]
        for k=1:ntimes #time loop
            for j=1:n^2 #individual loop
                ##############################
                #check death due to other causes



                ##############################
                #check neighbors states and sum of infected neighbors
                nn=Neigh[j,:]
                Ij=nn[findall(.!isnan.(nn))]
                nkj=sum(X[Ij]==2)
                M[1,1]=exp(-beta*nkj)
                M[2,1]=1-M[1,1]
                M[3,3]=sigma*M[1,1]
                M[2,3]=1-M[3,3]
                ##############################
                #update virus states, stacking probabilities method
                xi=rand
                pc=0
                jj=0
                while pc<xi
                    jj=jj+1
                    pc=pc+M[jj,X[j]]
                end
                Xnext[j]=jj
                ##############################
            end
            S[k,l]=sum(Xnext(X==1))
            Infected[k,l]=sum(Xnext(X==2))
            R[k,l]=sum(Xnext(X==3))
            X=Xnext
        end
end
#################################################################
#stochastic differential equation (SDE)
#arbirary num_s, num_r, c, d, P
#output tpoints mean, mean-std, mean+std
#start 11042019
#calculate propensities of reactions at time t (depends on N)
function vprop(N,c,d)
    v=zeros(num_r)
    for i=1:num_r
        index_n=1 #for mono-reaction
        index_find=false
        while index_find==false
            if d[i,index_n]<0
                index_find=true
                v[i]=c[i]*N[index_n]
            else
                index_n=index_n+1
            end
        end
    end
    return v
end

#calculate mean change
function change_mean(N,c,d)
    A=vprop(N,c,d)'*d
    return A
end
#calculate variance of change
function change_variance(N,c,d)
    B=zeros(Float16,num_s,num_s)
    for i=1:num_s
        for j=1:num_s
            sum=0
            for k=1:num_r
                sum=sum+vprop(N,c,d)[k]*d[k,i]*d[k,j]
            end
            B[i,j]=sum
        end
    end
    return B
end
#solving SDE using Euler Maruyama method
function StochDiffEquation(nMC,LEdims)
    start=time()
    pop_trace = zeros(Float64,num_s,num_t_SDE)
    MC_pop_trace = zeros(Float64,nMC, num_s,num_t_SDE)
    MC_pop_mean = zeros(Float64,num_s,num_t_SDE)
    MC_pop_std = zeros(Float64,num_s,num_t_SDE)
    MC_LE_mean = zeros(Float64,num_s,num_t_SDE)
    MC_LE_std = zeros(Float64,num_s,num_t_SDE)
    OS=zeros(Float64,nMC,num_t)
    OSmean=zeros(Float64,num_t)
    OSstd=zeros(Float64,num_t)
    LEsum=zeros(Float64,nMC)
    LEmean::Float64=0
    LEstd::Float64=0
    pop_trace[:,1] = pop_init
    dist_norm = Normal(0,1)
    for n=1:nMC
        for t=2:num_t_SDE
            dW=rand(dist_norm,num_s)
            A=change_mean(pop_trace[:,t-1],c,d)
            B=change_variance(pop_trace[:,t-1],c,d)
            Beigvals=eigvals(B)
            for m=1:num_s
                if Beigvals[m]<0
                    Beigvals[m]=0
                end
            end
            Bdiag=Diagonal(Beigvals)
            L=eigvecs(B)
            Linv=inv(L)
            Bsqrt=L*sqrt(Bdiag)*Linv
            pop_trace[:,t]=pop_trace[:,t-1]+A'*(dt/dt_scale)+Bsqrt*sqrt(dt/dt_scale)*dW
            #avoiding negative population
            #current approach (10/05/2019)
            for j=1:num_s
                if pop_trace[j,t]<0
                    pop_trace[j,t]=0
                end
            end
        end
        MC_pop_trace[n,:,:]=pop_trace

        for i=1:length(LEdims)
            LEsum[n]=LEsum[n]+sum(pop_trace[LEdims[i],:])
        end
        LEsum[n]=LEsum[n]/popsize
    end
    for i=1:num_s
        x=MC_pop_trace[:,i,:]
        MC_pop_mean[i,:]=mean(x;dims=1)
        MC_pop_std[i,:]=std(x;dims=1)
    end
    #return tpoints, MC_pop_mean, MC_pop_mean+MC_pop_std, MC_pop_mean-MC_pop_std
    #return tpoints, MC_pop_mean[:,index_pull], MC_pop_mean[:,index_pull]+MC_pop_std[:,index_pull], MC_pop_mean[:,index_pull]-MC_pop_std[:,index_pull]
    #calculate overall survival
    for j=1:nMC
        for i=1:length(LEdims)
            OS[j,:]=OS[j,:]+MC_pop_trace[j,LEdims[i],:]
        end
    end
    OSmean=mean(OS/popsize;dims=1)
    OSstd=std(OS/popsize;dims=1)
    #calculate life expectancy
    LEmean=mean(LEsum)
    LEstd=std(LEsum)
    #####
    elapsedtime=time()-start
    return elapsedtime, LEmean,  LEstd, tpoints, OSmean, OSstd, MC_pop_mean[:,index_pull], MC_pop_mean[:,index_pull]+MC_pop_std[:,index_pull], MC_pop_mean[:,index_pull]-MC_pop_std[:,index_pull]

end

# #################################################################
#microsimulation
#arbirary num_s, num_r, c, d, P
#output mean, mean-std, mean+std
#start 11042019
function MicroSimulation(nMC, LEdims)
    start=time()
    P=exp(Q)
    #check stochastic matrix, make sure row sum = 1
    last_index = zeros(Int8,num_s)
    row_sum = zeros(num_s)
    cum_sum = zeros(num_s,num_s)
    for i=1:num_s
        row_sum[i]=0
        for j=1:num_s
            if P[i,j]!==0 && row_sum[i] < 1
                row_sum[i]=row_sum[i]+P[i,j]
                cum_sum[i,j]=row_sum[i]
                last_index[i]=j
            else
                row_sum[i]=row_sum[i]
            end
        end
    end
    for i=1:num_s
        if sum(P[i,:]) !==1
            P[i,last_index[i]]=1-cum_sum[i,last_index[i]-1]
        end
    end
    print(P)
    @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # should be square
    # create vector of discrete RVs for each row
    dists = [Categorical(P[i, :]) for i in 1:N]
    # setup the simulation
    MC_pop_trace = zeros(Int32, nMC,popsize,num_t)
    pop_trace = zeros(Int32,popsize,num_t)
    MC_pop_count = zeros(Int32,nMC, num_s,num_t)
    MC_pop_mean = zeros(Float16,num_s,num_t)
    MC_pop_std = zeros(Float16,num_s,num_t)
    MC_LE_mean = zeros(Float64,num_s,num_t_SDE)
    MC_LE_std = zeros(Float64,num_s,num_t_SDE)
    LEsum=zeros(Float64,nMC)
    OS=zeros(Float64,nMC,num_t)
    OSmean=zeros(Float64,num_t)
    OSstd=zeros(Float64,num_t)
    LEmean::Float64=0
    LEstd::Float64=0
    pop_count = zeros(Int32, num_s,num_t)
    X = fill(0, num_t) # allocate memory, or zeros(Int64, sample_size)
    for n=1:nMC
        for i=1:popsize
            X[1] = 1 # set the pop_initial state
            for t in 2:num_t
                dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
                X[t] = rand(dist) # draw new value
            end
            pop_trace[i,:] = X
        end
        # collecting results: counting number of people in each state for each stoichiometries
        MC_pop_trace[n,:,:] = pop_trace
        for t=1:tfin
            for i=1:num_s
                pop_count[i,t] = count(pop_trace[:,t].== i)
                #pop_count[i,t] = count(pop_trace[pop_trace[:,t].== i,t])
            end
        end
        MC_pop_count[n,:,:]=pop_count

        for i=1:length(LEdims)
            LEsum[n]=LEsum[n]+sum(pop_count[LEdims[i],:])
        end
        LEsum[n]=LEsum[n]/(popsize*dt_scale)
    end
    for i=1:num_s
        x=MC_pop_count[:,i,:]
        MC_pop_mean[i,:]=mean(x;dims=1)
        MC_pop_std[i,:]=std(x;dims=1)
    end
    #calculate overall survival
    for j=1:nMC
        for i=1:length(LEdims)
            OS[j,:]=OS[j,:]+MC_pop_count[j,LEdims[i],:]
        end
    end
    OSmean=mean(OS/popsize;dims=1)
    OSstd=std(OS/popsize;dims=1)
    #calculate life expectancy
    LEmean=mean(LEsum)
    LEstd=std(LEsum)
    #####

    elapsedtime = time()-start
    return P, elapsedtime, LEmean, LEstd, tpoints, OSmean, OSstd, MC_pop_mean, MC_pop_mean+MC_pop_std, MC_pop_mean-MC_pop_std

end
#################################################################
#end of module
end
