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
#using DifferentialEquations, DiffEqBiological
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
function SetParamsBasic(X_names_in, Y_names_in, t0_in,tfin_in,dt_in,age_start_in, pop_init_in,popsize_in,rDeath_in)
     global X_names=X_names_in
     global Y_names=Y_names_in
     global t0 = t0_in #default is 0 (year)
     global tfin=tfin_in #year
     global dt=dt_in #virus timescale: default: 1-month
     global age_start = age_start_in #starting age
     global pop_init=pop_init_in
     global popsize=popsize_in
     global tpoints = collect(t0:dt:tfin)
     global num_t = length(tpoints)
     global index_pull=collect(1:dt*dt_scale:tfin*dt_scale+dt_scale-1)
     global rDeath=rDeath_in
end

#parameters virus
function SetParamsVirus(beta_in,sigma_in, rDeathVirusDis_in)
     global beta=beta_in #infectivity
     global sigma=sigma_in
     global rDeathVirusDis=rDeathVirusDis_in
end

#parameters cancer
function SetParamsCancer(rIncCancer_in, rIncCancerRR_in, rSoujournPre_in, rSoujourn_in, pStage_in, pStageScreen_in,
    rDeathCancerStage_in)
     global rIncCancer = rIncCancer_in #value must be less than 1
     global rIncCancerRR=rIncCancerRR_in
     global rIncCancerVirus = rIncCancer*rIncCancerRR
     global rSojournPre = rSojournPre_in
     global rSojourn = rSoujourn_in
     global pStage = pStage_in #make sure it's one-dimension
     global pStageScreen = pStageScreen_in #make sure it's one-dimension
     global rDeathCancerStage=rDeathCancerStage_in #vector of cancer death rates by stage (column)
end

function SetParamsScreen(pScreenSens_in, pScreenSpec_in, pScreenClin_in, intScreen)
     global pScreenSens = pScreenSens_in
     global pScreenSpec = pScreenSpec_in
     global pScreenClin = pScreenClin_in
     global intScreen = intScreen_in #in years
end
#################################################################
function simulate(nsims_in, screening_strategy_in)
#outputs dim: nsims x times
    nsims=nsims_in
    screening_strategy=screening_strategy_in
    n=sqrt(popsize)
    ntimes=num_t
    ##########################################################################################
    for l=1:nsims #monte carlo samples loop
        #time stamps for event for each individual
        #colums: 1: tDeath, 2: tIncVirus, 3: tIncCancer, 4: tDiagCancer,
        #5: tIncVirusDis, 6: tDeathCancer, 7:tDeathVirusDis
        tstamps = zeros(n^2,7)
        #current (clock) time
        tnow = t0
        age_now = age_start
        ################################
        #setting up neighboring states matrix
        I = reshape(1:n^2,n,n)
        Iext = [NaN for i=1:n+2 for j=1:n+2]
        Iext[2:n+1,2:n+1]=I
        Neigh=[NaN for i=1:n^2 for j=1:4]
        #west neighbors
        aux=Iext[2:n+1,1:n]
        Neigh[:,1]=aux'[:]
        #east neibhgors
        aux=Iext[2:n+1,3:n+2]
        Neigh[:,2]=aux'[:]
        #south neibhgors
        aux=Iext[3:n+2,2:n+1]
        Neigh[:,3]=aux'[:]
        #north neibhgors
        aux=Iext[1:n,2:n+1]
        Neigh[:,4]=aux'[:]
        #initial distribution of infected
        ninfect=10
        #state vectors for virus part
        X=ones[n^2,ntimes]
        X[1:2]=2; X[end]=2
        X[round(n^2/2)]=2
        X[round(n^2/3)]=2
        X[round(n^2/4)]=2
        X[round(n^2/5)]=2
        X[round(n^2/6)]=2
        X[round(n^2/7)]=2
        X[round(n^2/8)]=2
        Xnext=[NaN for i=1:n^2 for j=1:1]
        ################################
        #state trace for cancer and death parts: individuals*ntimes
        Y=ones[n^2,ntimes] #all start in no cancer
        #coding for state_names: 1: no cancer,
        #2: cancer-normal, 3: cancer-normal-diag-screen, 4: cancer-normal-diag-clin
        #5: cancer-virus, 6: cancer-virus-diag, 7: cancer-virus-diag-clin
        #8: death-other 9: death-virus 10: death-cancer
        ################################
        #flag time to cancer death prior to cancer incidence
        #draw time to cancer event point for all individuals (pre-allocated) - relative time
        #columns: time to cancer incidence-no virus, time to cancer incidence-infected, duration of
        distExpIncCancer=Exponential(rIncCancer)
        distExpSojourn=Exponential(rSojourn)
        distExpSojournPre=Exponential(rSojournPre)
        tsojourn=hcat(rand(distExpIncCancer,n^2),
        [0 for j=1:n^2]',
        [rand(distExpSojourn,n^2),
        [rand(distExpSojournPre,n^2))
        tsojourn[:,2]=rIncCancerRR*tsojourn[:,1] #make time to cancer due to virus sooner
        # tsojourn=hcat([Exponential(rIncCancer) for j=1:n^2]',
        # [0 for j=1:n^2]',
        # [Exponential(rSojourn) for j=1:n^2]',
        # [Exponential(rSojournPre) for j=1:n^2]')
        ################################
        #time stamps (time when event happens) - absolute time
        #colums: 1: tDeath, 2: tIncVirus, 3: tIncCancer, 4: tDiagCancer,
        #5: tIncVirusDis, 6: tDeathCancer, 7:tDeathVirusDis
        tstamps=hcat([0 for j=1:n^2]',
        [0 for j=1:n^2]',
        [0 for j=1:n^2]',
        [0 for j=1:n^2]',
        [0 for j=1:n^2]',
        [0 for j=1:n^2]',
        [0 for j=1:n^2]')
        ################################
        #draw time to death for all individuals (pre-allocated) - relative time
        #columns: death to OC, death to virus, death to cancer
        distExpDeathOC=Exponetial(rDeath)
        distExpDeathVirusDis=Exponential(rDeathVirusDis)
        tdeath=hcat(rand(distExpDeathOC,n^2),
        [9999 for j=1:n^2]',
        [9999 for j=1:n^2]')
        # tdeath=hcat([Exponential(rDeath) for j=1:n^2]',
        # [Exponential(rDeathVirusDis) for j=1:n^2]',
        # [9999 for j=1:n^2]')
        #flags for death (to execute certaint part of codes)
        deathflag=zeros(n^2,1)
        ################################
        #screening-related variables
        diagnosed=false #whether individual is diagnosed = removed from screening_strategy
        screening_schedule = collect(t0:intScreen*12:tfin)
        distMultinomStage = Multinomial(1,pStage) #set multinomial distribution for stage distribution at detection
        distMultinomStageScreen = Multinomial(1,pStageScreen)
        CancerStage = [NaN for mm=1:n^2] #storing cancer stage at detection during clinicaly-detectable phase
        CancerStageScreen = [NaN for mm=1:n^2] #storing cancer stage at detection during screen-detectable preclinical phase
        ##########################################################################################
        #start cohort simulation
        for k=1:ntimes #time loop (time progresses according to dt scale (default: month))
            #update time-related variables
            age_now=age_now+dt
            tnow=tnow+dt
            ##############################
            for j=1:n^2 #individual loop
                ##############################
                if deathflag==1 #flag to check whether j-th person is dead or not
                    Y[j,k+1]=Y[j,k]
                    X[j,k+1]=X[j,k]
                else #still alive
                    ############################################################
                    ###DEATH PART
                    ############################################################
                    #check death with hierarchy
                    if X[j,k]==2 #infected
                        if Y[j,k]==1 #no cancer
                            tdeath_min=argmin(tdeath[j,:])
                            if tnow<tdeath[j,tdeath_min]
                                sim_continue = true
                            else
                                Y[j,k+1]=7+tdeath_min #update state
                                sim_continue = false
                                deathflag[j]=1
                            end
                        elseif Y[j,k] > 1 & Y[j,k] <8 #with cancer
                            tdeath_min=argmin(tdeath[j,:])
                            if tnow<tdeath[j,tdeath_min]
                                sim_continue = true
                            else
                                Y[j,k+1]=7+tdeath_min #update state
                                sim_continue = false
                                deathflag[j]=1
                            end
                        end # cancer check
                    elseif X[j,k]==1 #no virus
                        if Y[j,k]==1 #no cancer
                            tdeath_min=argmin(tdeath[j,:])
                            if tnow<tdeath[j,tdeath_min]
                                sim_continue = true
                            else
                                Y[j,k+1]=7+tdeath_min #update state
                                sim_continue = false
                                deathflag[j]=1
                            end
                        elseif Y[j,k] > 1 & Y[j,k] <8 #with cancer
                            tdeath_min=argmin(tdeath[j,:]) #tdeath must include t deaht due to cancer (screening+stage dependent)
                            if tnow<tdeath[j,tdeath_min]
                                sim_continue = true
                            else
                                Y[j,k+1]=7+tdeath_min #update state
                                sim_continue = false
                                deathflag[j]=1
                            end
                        end # cancer check
                    else #person is already death, X=0
                        sim_continue = false
                    end #infectivity check
                    # if age_now > age_start+tstamps[j,1]
                    #     X[j]=0
                    #     tstamps[tstamps.=0].= NaN #flag other times if they are still "0"
                    # else
                    #check whether to continue (due to death)
                    if sim_continue==false #no need to update SIR or screening
                        Y[j,k+1]=Y[j,k]
                        X[j,k+1]=X[j,k]
                        Xnext[j]=X[j,k] # may be obsolete
                    else
                        ############################################################
                        ###CANCER DETECTION PART
                        ############################################################
                        if Y[j,k]==1 #no cancer
                            if X[j,k]==1 #no virus
                                if tnow > tsojourn[j,1] & tnow <tsojourn[j,1]+tsojourn[j,3] #cancer arises and in preclinical dp
                                    Y[j,k]=2 #is with cancer
                                    #no change of detection
                                elseif tnow > tsojourn[j,1]+tsojourn[j,3] & tnow <tsojourn[j,1]+tsojourn[j,3]+tsojourn[j,4] #cancer arises and in preclinical
                                    Y[j,k]=2 #is with cancer
                                ##############################
                                #check whether screening is on
                                    if screening_strategy==true
                                #check whether current time step is within screening interval
                                        if isempty(findall(x -> x.==tnow,screening_schedule))!=true
                                            #is detected?
                                            u=rand()
                                            if u<pScreenSens #cancer is detected by screening
                                                Y[k,j]=3
                                                #draw pStage
                                                stageindex=findall(x -> x.==1, vec(rand(distMultinomStage)))
                                                CancerStage[j]=stageindex #storing stage at diagnosis (screen)
                                                #time to death based on stage
                                                distExpDeathCancer=Exponential(stageindex)
                                                tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                            end
                                            #if not detected, Y remains at 2 (has cancer but undetected)
                                        end #end screeening interval checking
                                    end #screening_strategy condition
                                else #cancer in clinically-detectable phase
                                    Y[j,k]=2 #is with cancer
                                    u=rand()
                                    if u<pScreenClin #cancer is clinically detected
                                        Y[k,j]=4
                                        #draw pStage
                                        stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                        CancerStageScreen[j]=stageindex #storing stage at diagnosis (clinical)
                                        #time to death based on stage
                                        distExpDeathCancer=Exponential(stageindex)
                                        tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                    end
                                end #end sojourn preclinical ND
                            elseif X[j,k]==2 #with virus
                                if tnow > tsojourn[j,2] & tnow <tsojourn[j,2]+tsojourn[j,3]#cancer arises and in preclinical dp
                                    Y[j,k]=5 #is with cancer
                                    #no change of detection
                                elseif tnow > tsojourn[j,2]+tsojourn[j,3] & tnow <tsojourn[j,2]+tsojourn[j,3]+tsojourn[j,4] #cancer arises and in preclinical
                                    Y[j,k]=5 #is with cancer
                                ##############################
                                    if screening_strategy==true
                                #check whether current time step is within screening interval
                                        if isempty(findall(x -> x.==tnow,screening_schedule))!=true
                                            #is detected?
                                            u=rand()
                                            if u<pScreenSens #cancer is detected by screening
                                                Y[k,j]=6
                                                #draw pStage
                                                stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                                CancerStageScreen[j]=stageindex
                                                #time to death based on stage
                                                distExpDeathCancer=Exponential(stageindex)
                                                tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                            end
                                        end #end screeening interval checking
                                    end #screening strategy condition
                                else #cancer in clinically-detectable phase
                                    Y[j,k]=5 #is with cancer
                                    u=rand()
                                    if u<pScreenClin #cancer is clinically detected
                                        Y[k,j]=7
                                        #draw pStage
                                        stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                        CancerStageScreen[j]=stageindex #storing stage at diagnosis (clinical)
                                        #time to death based on stage
                                        distExpDeathCancer=Exponential(stageindex)
                                        tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                    end
                                end #end sojourn preclinical ND
                            end #end virus condition
                        elseif Y[j,k]==2 | Y[j,k]==5 #has cancer and not detected
                            if X[j,k]==1 #no virus
                                if tnow > tsojourn[j,1]+tsojourn[j,3] & tnow <tsojourn[j,1]+tsojourn[j,3]+tsojourn[j,4] #cancer arises and in preclinical
                                    #Y[j,k]=2 #is with cancer
                                ##############################
                                    if screening_strategy==true
                                    #check whether current time step is within screening interval
                                        if isempty(findall(x -> x.==tnow,screening_schedule))!=true
                                            #is detected?
                                            u=rand()
                                            if u<pScreenSens #cancer is detected by screening
                                                Y[k,j]=3
                                                #draw pStage
                                                stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                                CancerStageScreen[j]=stageindex
                                                #time to death based on stage
                                                distExpDeathCancer=Exponential(stageindex)
                                                tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                            end
                                        end #end screeening interval checking
                                    end #screening strategy condition
                                else #cancer in clinically-detectable phase
                                    #Y[j,k]=2 #is with cancer
                                    u=rand()
                                    if u<pScreenClin #cancer is clinically detected
                                        Y[k,j]=4
                                        #draw pStage
                                        stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                        CancerStageScreen[j]=stageindex #storing stage at diagnosis (clinical)
                                        #time to death based on stage
                                        distExpDeathCancer=Exponential(stageindex)
                                        tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                    end
                                end #end sojourn preclinical ND
                            elseif X[j,k]==2 #with virus
                                if tnow > tsojourn[j,2]+tsojourn[j,3] & tnow <tsojourn[j,2]+tsojourn[j,3]+tsojourn[j,4] #cancer arises and in preclinical
                                    #Y[j,k]=5 #is with cancer
                                ##############################
                                    if screening_strategy==true
                                    #check whether current time step is within screening interval
                                        if isempty(findall(x -> x.==tnow,screening_schedule))!=true
                                            #is detected?
                                            u=rand()
                                            if u<pScreenSens #cancer is detected by screening
                                                Y[k,j]=6
                                                #draw pStage
                                                stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                                CancerStageScreen[j]=stageindex
                                                #time to death based on stage
                                                distExpDeathCancer=Exponential(stageindex)
                                                tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                            end
                                        end #end screeening interval checking
                                    end #screening strategy condition
                                else #cancer in clinically-detectable
                                    #Y[j,k]=5 #is with cancer
                                    u=rand()
                                    if u<pScreenClin #cancer is clinically detected
                                        Y[k,j]=7
                                        #draw pStage
                                        stageindex=findall(x -> x.==1, vec(rand(distMultinomStageScreen)))
                                        CancerStageScreen[j]=stageindex #storing stage at diagnosis (clinical)
                                        #time to death based on stage
                                        distExpDeathCancer=Exponential(stageindex)
                                        tdeath[j,3]=tnow+rand(distExpDeathCancer) #store time to  death for used in above death hierarchy section
                                    end
                                end #end sojourn preclinical ND
                            end #end virus condition
                        end #cancer incidence condition
                        ############################################################
                        ###VIRUS PART
                        ############################################################
                        #check neighbors states and sum of infected neighbors
                        nn=Neigh[j,:]
                        Ij=nn[findall(.!isnan.(nn))]
                        nkj=sum(X[Ij]==2) # 2 is infected #0 is death
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
                        Xnext[j]=jj #may be obsolete
                        X[j,k+1]=jj
                        Y[j,k+1]=Y[j,k]
                        ###############
                        #draw death due to virus, if status changed to infected for the 1st time
                        if X[j,k]==1 & X[j,k+1]==2
                            tdeath[j,2]=tnow+rand(distExpDeathVirusDis)
                        end
                    end #end of sim_continue condition
                end #end of deathflag condition
            end #end of person loop
            #####################################
            S[k,l]=sum(Xnext(X[:,l]==1))
            Infected[k,l]=sum(Xnext(X==2))
            R[k,l]=sum(Xnext(X==3))
            #X[:,]=Xnext
            #####################################
            #update neighboring states based on updated X
            I = reshape(X,n,n)
            Iext = [NaN for i=1:n+2 for j=1:n+2]
            Iext[2:n+1,2:n+1]=I'
            Neigh=[NaN for i=1:n^2 for j=1:4]
            #west neighbors
            aux=Iext[2:n+1,1:n]
            Neigh[:,1]=aux'[:]
            #east neibhgors
            aux=Iext[2:n+1,3:n+2]
            Neigh[:,2]=aux'[:]
            #south neibhgors
            aux=Iext[3:n+2,2:n+1]
            Neigh[:,3]=aux'[:]
            #north neibhgors
            aux=Iext[1:n,2:n+1]
            Neigh[:,4]=aux'[:]
        end #end of time loop
    end #end of MC loop
end #end of simulate function

#################################################################

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
