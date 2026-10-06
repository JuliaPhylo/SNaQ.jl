## test to see if we get lik=0.0 with perfect data
## with a network of 15 taxa, 3 hyb
## Claudia September 2016

net = readnewick("(15,(1,((14,(#H1,(((12,13),(11,#H3)),(7,((10)#H3,(8,9)))))),((((2,3))#H2,(6,(5,(#H2,4)))))#H1)));");
for E in net.edge E.length = 0.25 end
for H in net.hybrid getparentedge(H).gamma = 0.35; getparentedgeminor(H).gamma = 0.65; end
df = computeexpectedDataCF(net)

val = SNaQ.optimizetopology!(net, df)

@test val ≈ 0.0 atol=1e-10 # || error("not correct likelihood with perfect data")
