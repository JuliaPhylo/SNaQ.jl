using PhyloNetworks, Test, StatsBase, Random
import SNaQ:
    move_reticulate_target!, move_reticulate_origin!,
    semidirect_network!

#### Helper functions ####
function reload_labelled_net()
    net = readnewick("((a,b)i1,((c,#H1)i6,(d,((e,f)i5)#H1)i3)i2);")
    semidirect_network!(net)
    net.hybrid[1].name = "i4"   # do it here so that PhyloNetworks doesn't through a warning
    return net
end
get_nodes(net::HybridNetwork, names::String...) = [
    net.node[findfirst(node -> node.name == name, net.node)] for name in names
]




#############################################################################
# FIRST: THESE CASES HAVE CAUSED BUGS IN THE PAST - MAKE SURE THEY WORK NOW #
#############################################################################

net = readnewick("(((c1:0.5,c2:0.5)i3:0.622,((e:0.5,f:0.5)i10:0.206)#Hi4:0.206::0.739)i4:0.31,((d2:0.25,(((#Hi4:0.206::0.261,d1:0.5)i7:0.008)#H2:0.008::0.978)#H3:0.008::0.774)i5440:0.662,#H3:0.008::0.226)i1976:0.662,(b:0.5,(a:0.25,#H2:0.008::0.022)i12:0.709)i15:1.004)i16;");
net.isrooted = false;
H, i7, Hi4 = get_nodes(net, "H3", "i7", "Hi4");
E = [e for e in i7.edge if Hi4 in e.node][1];
move_reticulate_target!(net, H, E)
@test all([getparentedge(hyb).gamma + getparentedgeminor(hyb).gamma ≈ 1 for hyb in net.hybrid])


net = readnewick("(((b:0.75)#H3:0.5::0.972,a:0.96)i15:1.081,(((e:0.5,f:0.5)i10:0.208)#Hi4:0.208::0.72,((d1:0.5,d2:0.75)i7:0.095,#H4:0.5::0.007)i1817:0.549)i12:0.281,((#Hi4:0.208::0.28,(c2:0.447)#H4:0.5::0.993)i3:0.0,(c1:0.132,#H3:0.5::0.028)i5645:0.993)i4:1.127)i16;");
net.isrooted = false;
H, i12, i1817 = get_nodes(net, "Hi4", "i12", "i1817");
E = [e for e in i12.edge if i1817 in e.node][1];
move_reticulate_origin!(net, H, E)
@test all([getparentedge(hyb).gamma + getparentedgeminor(hyb).gamma ≈ 1 for hyb in net.hybrid])


net = readnewick("((c2:1.162,(((a:1.136,(b:0.548,#H5:0.5::0.043)i8475:0.625)H3:1.187,((d1:1.038,d2:1.301)i10:0.35,#H4:0.5::0.32)i9744:0.539)i3628:0.149,((f:0.776,(#H6:0.5::0.465)#H4:0.5::0.68)i5:0.481,(e:0.15)#H6:0.5::0.535)i3056:0.481)i12:0.831)i2450:0.5,(c1:0.285)#H5:0.5::0.957)i7;");
net.isrooted = false;
H4, i12, i3056 = get_nodes(net, "H4", "i12", "i3056");
E = [e for e in i12.edge if i3056 in e.node][1];
move_reticulate_target!(net, H4, E)
@test all([getparentedge(hyb).gamma + getparentedgeminor(hyb).gamma ≈ 1 for hyb in net.hybrid])


