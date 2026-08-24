# Occasionally edge lengths of imported networks will be -1 (e.g. edges above leaves,
# or edges above hybrid nodes that were inferred with previous versions of SNaQ).
#
# Previously there was a bug where in the latter case an error would be thrown for
# negative edge lengths, because these previously unidentifiable edges are now identifiable.
# This test file checks:
# 1. Hybrid edges with negative length are automatically changed to have length 0.0 when
#	utilized in functions like `computeSNaQscore!` or `fitnumericalparameters!`
# 2. Leaf edges are allowed to have negative length without issue (because they are
#	still unidentifiable)


global dcf0::DataCF;
global dcf1::DataCF;
global dcf2::DataCF;

# set edges to have different, fixed values so that we can compare the
# equality of the generate DataCF expected CFs
sett!(edge::PhyloNetworks.Edge) = (edge.length = edge.hybrid ? 0.0 : Float64(abs(edge.number / 10.0)))
net = readnewick("(((a,#H2),(b,#H1)),(((c,d),((e,(f)#H1))#H2)));");

# Fully specified edge lengths
for E in net.edge
	E.length = sett!(E)
	E.gamma = E.hybrid ? 0.5 : -1
end
@test (global dcf0 = computeexpectedDataCF(net); true)	# passes as long as there are no errors thrown


# Unspecified leaf edge lengths
for E in net.edge
	E.length = getchild(E).leaf ? -1 : sett!(E)
end
@test (global dcf1 = computeexpectedDataCF(net); true)


# Unspecified leaf and hybrid edge lengths
for E in net.edge
	E.length = getchild(E).leaf || getchild(E).hybrid ? -1 : sett!(E)
end
@test (global dcf2 = computeexpectedDataCF(net); true)


# All computed expected CFs are identical
@test all(i -> all(dcf0.quartet[i].expCF .== dcf1.quartet[i].expCF .== dcf2.quartet[i].expCF), eachindex(dcf0.quartet))


# An error should still be thrown if there is a negative edge length that is not allowed
net.edge[10].length = -1.0
@test_throws ErrorException computeexpectedDataCF(net)
@test length(getnegativeedges(net)) == 1 && getnegativeedges(net)[1].number == 10
