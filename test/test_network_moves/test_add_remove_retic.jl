using PhyloNetworks, Test


################################################
# FIRST: WE KNOW THESE MOVES SHOULD BE INVALID #
################################################

net = readnewick("(a,((b,#H1)i2,(c,((d,e)i5)#H1)i4)i3)i1;")
e2 = getparentedge(net.leaf[4])
e1 = getparentedge(getparent(e2))
@test is_valid_add_hybrid(e1, e2, net)
@test !is_valid_add_hybrid(e2, e1, net)


################################################################
# SECOND: add_hybrid! SHOULD BE REVERSIBLE WITH remove_hybrid! #
################################################################

net = readnewick("(a,((b,#H1)i2,(c,((e,d)i5)#H1)i4)i3)i1;");
net0 = deepcopy(net);

e1 = getparentedge(net.leaf[4])
e2 = getparentedge(getparent(e1))

newH = add_hybrid!(e2, e1, net);
remove_hybrid!(newH, net)

# (e,d) gets swapped to (d,e) after the operation, but the Newick is equivalent otherwise
@test replace(writenewick(net), "d" => "e") == replace(writenewick(net0), "d" => "e")






net = readnewick("(((d1:1.602,d2:1.602)i10:0.823,(b:0.463,(a:0.573)#H1:0.5::0.973)i3442:1.191)i13:0.111,((c2:0.757,c1:1.554)i6:0.799,(e:0.776,f:0.776)i5:0.497)i12:0.0,#H1:0.5::0.027)i7;"); net.isrooted = false;
@test_nowarn move_random_reticulate_origin!(net);

new_origin = getparentedge(net.leaf[3])
H = net.hybrid[1]
@test is_valid_move_reticulate_origin(H, new_origin, net)


net = readnewick("(b:1.095,(a:0.268,#H2:0.562::0.5)i3649:0.5,((((e:0.679,f:0.515)i5:0.496,(c2:0.51,c1:0.0)i13:0.8)i6:0.109,(d1:0.65,d2:1.602)H1:0.825)i10737:0.562)#H2:0.562::0.5)i7;");
net.isrooted = false;

try
    move_random_reticulate_origin!(net)
    @test false
catch e    buf = IOBuffer()
    showerror(buf, e)
    message = String(take!(buf))
    @test message == "No valid `move_random_reticulate_origin!` moves."
end
