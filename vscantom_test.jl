######################################################
#file to test CohMod module
include("./vscantom.jl")
using Plots
using GraphPlot
using LightGraphs
#testing vscantom.jl module
#Rowan Iskandar (rowan.iskandar@gmail.com)
#start 01092019
#two examples
######################################################
#coding for X state_names: 0: death, 1: susceptible, 2: virus
X_names_in=["death","S","V"]
#coding for Y state_names: 1: no cancer,
#2: cancer-normal, 3: cancer-normal-diag-screen, 4: cancer-normal-diag-clin
#5: cancer-virus, 6: cancer-virus-diag, 7: cancer-virus-diag-clin
#8: death-other 9: death-virus 10: death-cancer
Y_names_in=["no cancer", "cancer-noV-noD","cancer-noV-D-S", "cancer-noV-D-C",
"cancer-V-noD","cancer-V-D-S", "cancer-V-D-C",
"death-OC","death-V", "death-C"]
 t0_in = 0 #default is 0 (year)
 tfin_in=70 #year
 dt_in=1 #virus timescale: default: 1-month
 age_start_in = 30 #starting age
 pop_init_in=10000
 popsize_in=10000
 rDeath_in=0.001
vscantom.SetParamsBasic(X_names_in, Y_names_in,t0_in,tfin_in,dt_in,age_start_in, pop_init_in,popsize_in,rDeath_in)

#parameters virus
beta_in=1 #infectivity
sigma_in=0
rDeathVirusDis_in=0.01
SetParamsVirus(beta_in,sigma_in,rDeathVirusDis_in)



#parameters cancer
 rIncCancer_in = 0.7 #value must be less than 1
 rIncCancerRR_in=0.01
 rIncCancerVirus_in = rIncCancer*rIncCancerRR
 rSojournPre_in = 0.02
 rSojourn_in = 0.02
 pStage_in = [0.7 0.2 0.095 0.005] #make sure it's one-dimension
 pStageScreen_in = [0.8 0.15 0.049 0.001] #make sure it's one-dimension
 rDeathCancerStage_in=[0.001 0.005 0.01 0.1] #vector of cancer death rates by stage (column)

SetParamsCancer(rIncCancer_in, rIncCancerRR_in, rSoujournPre_in, rSoujourn_in, pStage_in, pStageScreen_in,
    rDeathCancerStage_in)


pScreenSens_in = 0.9
pScreenSpec_in = 0.97
pScreenClin_in = pScreenClin_in
intScreen_in = 1 #in years

SetParamsScreen(pScreenSens_in, pScreenSpec_in, pScreenClin_in, intScreen)

#################################################################
virus_prevention_in = false
screening_strategy_in = true
nsims_in = 1
simulate(nsims_in, virus_prevention_in, screening_strategy_in)
