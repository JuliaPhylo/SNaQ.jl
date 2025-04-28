# test function from a debugging case in simScript.jl
# in julia/debug/bigSimulation/n6/julia
# seed 4545
# Claudia August 2015

# to debug problem
tree="(3,(2,(((6,(5)#H9:0.91507):0.93066,(4,#H9:0.0):0.73688):0.0)#H7:1.79104::0.99498):0.11675,(1,#H7:0.04487::0.00502):0.4897);"
net0=readnewicklevel1(tree);
#printedges(net0)
gammaz!(net0.node[6], 1.0)  #-5
gammaz!(net0.node[8], 0.067) #-7

q1 = Quartet(1,["6","5","4","1"],[0.5,0.4,0.1]);
qnet = extractQuartet!(net0,q1);
#printedges(qnet)

identifyQuartet!(qnet)
qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
#printedges(qnet)

eliminateHybridization!(qnet)
qnet.which != 2 ? error("qnet which not correctly assigned") : nothing
qnet.numhybrids == 1 || error("qnet not correctly eliminated hyb")
k(qnet.hybrid[1]) == 4 || error("qnet not correclty identified")
typeHyb(qnet.hybrid[1]) == 5 || error("qnet not correclty identified")
isBadDiamondI(qnet.hybrid[1]) || error("qnet forgot it is a bad diamondI")

updateSplit!(qnet)
qnet.split != [-1,-1,-1,-1] ? error("qnet.split not correctly assigned") : nothing

updateFormula!(qnet)
qnet.formula != [-1,-1,-1] ? error("qnet.formula not correctly assigned") : nothing

calculateExpCF!(qnet)
qnet.expCF[2] < 0. || error("expCF not correctly calculated")

