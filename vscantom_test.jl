######################################################
#file to test CohMod module
include("./vscantom.jl")
using Plots
using GraphPlot
using LightGraphs
#testing CohMod.jl module
#Rowan Iskandar (rowan.iskandar@gmail.com)
#start 10042019
#two examples
######################################################
#example 1: 4-state example in
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205543
#set up parameters
c12 = 0.05
c13 = 0.01
c14 = 0.001
c23 = 0.1
c24 = 0.05
c34 = 1
popsize = 1000
nMC = 10000 # number of micro sim iterations
t0 = 0 #initial time
tfin = 50 #final time
dt = 1 #time step
dt_scale = 1
c = [c12 ;c13 ;c14 ;c23 ;c24 ;c34 ]
pop_init = [popsize 0 0 0] #initial state configuration
#generator matrix

state_names=["S1" "S2" "S3" "S4"]
#call module and set global parameter values
CohMod.SetParams(state_names,t0,tfin,dt,dt_scale,Q,c,d,pop_init,popsize)

@time begin
tpoints,  SDEmean, SDEnegsd, SDEpossd = CohMod.StochDiffEquation(nMC)
end
plot(tpoints,SDEmean[1,:])
plot!(tpoints,SDEmean[2,:])
plot!(tpoints,SDEmean[3,:])
plot!(tpoints,SDEmean[4,:])
#microsimulation
@time begin
tpoints, MCmean, MCnegsd, MCpossd  = CohMod.MicroSimulation(nMC)
end
plot(tpoints,MCmean[1,:])
plot!(tpoints,MCmean[2,:])
plot!(tpoints,MCmean[3,:])
plot!(tpoints,MCmean[4,:])
#ODE
@time begin
ODE= CohMod.DiffEquation( true,false,"accurate")
end
plot(ODE)
#CMEgil
@time begin
MEgil_out = CohMod.MasterEqGil(1)
end
