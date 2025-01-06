# functions to read/write level-1 networks

# function to clean topology after readTopology
# looks for:
# TREE:
# - all tree edges must have gamma=1. fixit: cannot point out which doesn't,
#   only shows error.
# - internal nodes with only 2 edges and solves this case.
# - polytomies and choose one resolution at random, issuing a warning
# NETWORK:
# - number of hybrid edges per hybrid node:
#   if 0,1: error (with warning in old functions)
#   if >2: error of hybrid polytomy
#   if 2: check number of tree edges
# - number of tree edges per hybrid node:
#   if 0: leaf hybrid, add child
#   if >1: expand child
#   if 1: check values of gamma:
# - gammas: need to sum to one and be present.
#   error if they do not sum up to one
#   default values of 0.1,0.9 if not present
# leaveRoot=true: leaves the root even if it has only 2 edges (for plotting), default=false
function cleanAfterRead!(net::HybridNetwork, leaveRoot::Bool)
    # set e.containroot to !e.hybrid: updated later by updateAllReadTopology as required by snaq!
    for e in net.edge e.containroot = !e.hybrid; end
    nodes = copy(net.node)
    for n in nodes
        if isNodeNumIn(n,net.node) # very important to check
            if size(n.edge,1) == 2 # delete n if:
                if (!n.hybrid && (!leaveRoot || !isEqual(getroot(net),n)) ||
                    (n.hybrid && sum(e.hybrid for e in n.edge) == 1))
                    deleteIntNode!(net,n)
                    continue # n was deleted: skip the rest
                end
            end
            if !n.hybrid
                if size(n.edge,1) > 3
                    @debug "warning: polytomy found in node $(n.number), random resolution chosen"
                    resolvetreepolytomy!(net,n);
                end
                hyb = count([e.hybrid for e in n.edge])
                if hyb == 1
                    hasHybEdge(n) == true;
                elseif hyb > 1
                    @warn "strange tree node $(n.number) with more than one hybrid edge, intersecting cycles maybe"
                end
            else
                hyb = count([e.hybrid for e in n.edge]);
                tre = length(n.edge) - hyb
                if hyb > 2
                    error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
                elseif hyb == 1
                    hybnodes = count([n.hybrid for n in net.node]);
                    if hybnodes == 1
                        error("only one hybrid node number $(n.number) with name $(net.names[n.number]) found with one hybrid edge attached")
                    else
                        error("current hybrid node $(n.number) with name $(net.names[n.number]) has only one hybrid edge attached. there are other $(hybnodes-1) hybrids out there but this one remained unmatched")
                    end
                elseif hyb == 0
                    @warn "hybrid node $(n.number) is not connected to any hybrid edges, it was transformed to tree node"
                    n.hybrid = false;
                else # 2 hybrid edges
                    if tre == 0 #hybrid leaf
                        @warn "hybrid node $(n.number) is a leaf, so we add an extra child"
                        addChild!(net,n);
                    elseif tre > 1
                        @warn "hybrid node $(n.number) has more than one child so we need to expand with another node"
                        expandChild!(net,n);
                    end
                    suma = sum([e.hybrid ? e.gamma : 0.0 for e in n.edge]);
                    # synchronizepartnersdata! already made suma ≈ 1.0, when non-missing,
                    # and already synchronized ismajor, even when γ's ≈ 0.5
                    if suma == -2.0 # hybrid edges have no gammas in newick description
                        println("hybrid edges for hybrid node $(n.number) have missing gamma's, set default: 0.9,0.1")
                        for e in n.edge
                            if e.hybrid
                                e.gamma = (e.ismajor ? 0.9 : 0.1)
                            end
                        end
                    end
                end
            end
        end
    end
    for e in net.edge
        if e.hybrid
          if e.gamma < 0.0 || e.gamma > 1.0 # no -1.0 (missing) γ's by now
            error("hybrid edge $(e.number) with γ = $(e.gamma) outside of [0,1]")
          end
        else
          e.gamma == 1.0 ||
            error("tree edge $(e.number) with γ = $(e.gamma) instead of 1.0")
        end
    end
end

cleanAfterRead!(net::HybridNetwork) = cleanAfterRead!(net,false)

# function to update the read topology after reading
# it will go over the net.hybrid array and check each
# of the hybridization events defined to update:
# - in cycle
# - contain root
# - gammaz
# it uses updateAllNewHybrid! function that
# returns: success (bool), hybrid, flag, nocycle, flag2, flag3
# if tree read, also check that contain root is true for all, ishybrid and hashybedge is false
# warning: needs to have run storeHybrids! before
# warning: it will stop when finding one conflicting hybrid
function updateAllReadTopology!(net::HybridNetwork)
    if(isTree(net))
        #@warn "not a network read, but a tree as it does not have hybrid nodes"
        all((e->e.containroot), net.edge) ? nothing : error("some tree edge has contain root as false")
        all((e->!e.hybrid), net.edge) ? nothing : error("some edge is hybrid and should be all tree edges in a tree")
        all((n->!hasHybEdge(n)), net.node) ? nothing : error("some tree node has hybrid edge true, but it is a tree, there are no hybrid edges")
    else
        if(!cleaned(net))
            for n in net.hybrid
                success,hyb,flag,nocycle,flag2,flag3 = updateAllNewHybrid!(n,net,false,true,false)
                if(!success)
                    @warn "current hybrid $(n.number) conflicts with previous hybrid by intersecting cycles: $(!flag), nonidentifiable topology: $(!flag2), empty space for contain root: $(!flag3), or does not create a cycle (probably problem with the root placement): $(nocycle)."
                    #cleaned!(net, false)
                end
            end
            @debug "before update partition"
            @debug begin printpartitions(net); "printed partitions" end
            for n in net.hybrid #need to updatePartition after all inCycle
                nocycle, edgesInCycle, nodesInCycle = identifyInCycle(net,n);
                updatePartition!(net,nodesInCycle)
                @debug begin printpartitions(net);
                    "partitions after updating partition for hybrid node $(n.number)"
                end
            end
        end
    end
end

# cleanAfterReadAll includes all the step to clean a network after read
function cleanAfterReadAll!(net::HybridNetwork, leaveRoot::Bool)
    @debug "check for 2 hybrid edges at each hybrid node -----"
    check2HybEdges(net)
    @debug "cleanBL -----"
    cleanBL!(net)
    @debug "cleanAfterRead -----"
    cleanAfterRead!(net,leaveRoot)
    @debug "updateAllReadTopology -----"
    updateAllReadTopology!(net) #fixit: it could break if leaveRoot = true (have not checked it), but we need to updateContainRoot
    if(!leaveRoot)
        @debug "parameters -----"
        parameters!(net)
    end
    @debug "check root placement -----"
    checkRootPlace!(net)
    getroot(net).leaf && @warn "root node $(getroot(net).number) is a leaf, so when plotting net, it can look weird"
    cleaned!(net, true) #fixit: set to false inside updateAllReadTopology if problem encountered
    net.isrooted = false
end

cleanAfterReadAll!(net::HybridNetwork) = cleanAfterReadAll!(net,false)

# function to read a topology from file name/tree directly and update it
# by calling updateAllReadTopology after
# leaveRoot=true if the root will not be deleted even if it has only 2 edges
# used for plotting (default=false)
# warning: if leaveRoot=true, net should not be used outside plotting, things will crash
function readTopologyUpdate(file::AbstractString, leaveRoot::Bool,verbose::Bool)
    @debug "readnewick -----"
    net = readnewick(file,verbose)
    cleanAfterReadAll!(net,leaveRoot)
    return net
end
readTopologyUpdate(file::AbstractString) = readTopologyUpdate(file, false, true)
readTopologyUpdate(file::AbstractString,verbose::Bool) = readTopologyUpdate(file, false, verbose)


"""
    readnewicklevel1(filename)
    readnewicklevel1(parenthetical format)

Similarly to [`PhyloNetworks.readnewick`](): read a tree or network in parenthetical
format, but this function enforces the necessary conditions for any
starting topology in SNaQ: non-intersecting cycles, no polytomies,
unrooted. It sets any missing branch length to 1.0,
and reduces any branch length above 10 to 10 (as it assumes branch lengths are in coalescent units).

If the network has a bad diamond II (in which edge lengths are γ's are not identifiable)
and if the edge below this diamond has a length `t` different from 0, then this length is
set back to 0 and the major parent hybrid edge is lengthened by `t`.
"""
readnewicklevel1(file::AbstractString) = readTopologyUpdate(file, false, true)

# to read multiple topologies: readmultinewicklevel1 is defined in readquartetdata.jl
# It calls readTopologyUpdate defined here, for level 1 networks.

# aux function to send an error if the number of hybrid attached to every
# hybrid node is >2
function check2HybEdges(net::HybridNetwork)
    for n in net.hybrid
        hyb = sum([e.hybrid for e in n.edge]);
        if hyb > 2
            error("hybrid node $(n.number) has more than two hybrid edges attached to it: polytomy that cannot be resolved without intersecting cycles.")
        end
    end
end


# aux function to check if the root is placed correctly, and re root if not
# warning: it needs updateContainRoot set
function checkRootPlace!(net::HybridNetwork; verbose=false::Bool, outgroup="none"::AbstractString)
    if(outgroup == "none")
        if !canBeRoot(getroot(net))
            verbose && println("root node $(getroot(net).number) placement is not ok, we will change it to the first found node that agrees with the direction of the hybrid edges")
            for i in 1:length(net.node)
                if(canBeRoot(net.node[i]))
                    net.rooti = i
                    break
                end
            end
        end
    else # outgroup
        tmp = findall(n -> n.name == outgroup, net.leaf)
        if length(tmp)==0
            error("leaf named $(outgroup) was not found in the network.")
        elseif length(tmp)>1
            error("several leaves were found with name $(outgroup).")
        end
        leaf = net.leaf[tmp[1]]
        leaf.leaf || error("found outgroup not a leaf: $(leaf.number), $(outgroup)")
        length(leaf.edge) == 1 || error("found leaf with more than 1 edge: $(leaf.number)")
        other = getOtherNode(leaf.edge[1],leaf);
        if(canBeRoot(other))
            net.rooti = getIndexNode(other.number,net)
        else
            throw(RootMismatch("outgroup $(outgroup) contradicts direction of hybrid edges"))
        end
    end
    canBeRoot(getroot(net)) || error("tried to place root, but couldn't. root is node $(getroot(net))")
end



# see full docstring below
# Need HybridNetwork input, since QuartetNetwork does not have root.
function writenewick_level1(net0::HybridNetwork, di::Bool, str::Bool, namelabel::Bool,outgroup::AbstractString, printID::Bool, roundBL::Bool, digits::Integer, multall::Bool)
    s = IOBuffer()
    writenewick_level1(net0,s,di,namelabel,outgroup,printID,roundBL,digits, multall)
    if str
        return String(take!(s))
    else
        return s
    end
end

# warning: I do not want writenewick_level1 to modify the network if outgroup is given! thus, we have updateRoot, and undoRoot
# note that if printID is true, the function is modifying the network
function writenewick_level1(net0::HybridNetwork, s::IO, di::Bool, namelabel::Bool,
           outgroup::AbstractString, printID::Bool, roundBL::Bool, digits::Integer, multall::Bool)
    global CHECKNET
    net = deepcopy(net0) #writenewick_level1 needs containroot, but should not alter net0
    # numBad(net) == 0 || println("network with $(numBad(net)) bad diamond I. Some γ and edge lengths t are not identifiable, although their γ * (1-exp(-t)) are.")
    if printID
        setNonIdBL!(net) # changes non identifiable BL to -1.0, except those in/below bad diamonds/triangles.
    end
    assignhybridnames!(net)
    if net.numnodes == 1
        print(s,string(getroot(net).number,";")) # error if 'string' is an argument name.
    else
        if(!isTree(net) && !cleaned(net))
            @debug "net not cleaned inside writenewick_level1, need to run updateContainRoot"
            for n in net.hybrid
                flag,edges = updateContainRoot!(net,n)
                flag || error("hybrid node $(n.hybrid) has conflicting containroot")
            end
        end
        updateRoot!(net,outgroup)
        #@debug begin printEverything(net); "printed everything" end
        CHECKNET && canBeRoot(getroot(net))
        if(multall)
            mergeLeaves!(net)
            ## make sure the root is not on a leaf
            ## This is a band aid: need to check the order of write/root/merge leaves on multiple allele cases
            if outgroup != "none"
                try
                    checkRootPlace!(net,outgroup=outgroup) ## keeps all attributes
                catch err
                    if isa(err, RootMismatch)
                        println("RootMismatch: ", err.msg,
                                """\nThe estimated network has hybrid edges that are incompatible with the desired outgroup.
                        Reverting to an admissible root position.
                        """)
                    else
                        println("error trying to reroot: ", err.msg);
                    end
                    checkRootPlace!(net,verbose=false) # message about problem already printed above
                end
            else
                checkRootPlace!(net,verbose=false) #leave root in good place after snaq
            end
        end
        writesubtree!(s, net, di,namelabel, roundBL,digits,true)
    end
    # outgroup != "none" && undoRoot!(net) # not needed because net is deepcopy of net0
    # to delete 2-degree node, for snaq.
end

writenewick_level1(net::HybridNetwork,di::Bool,str::Bool,namelabel::Bool,outgroup::AbstractString,printID::Bool) = writenewick_level1(net,di,str,namelabel,outgroup,printID, false,3, false)
# above: default roundBL=false (at unused digits=3 decimal places)
writenewick_level1(net::HybridNetwork,printID::Bool) = writenewick_level1(net,false, true,true,"none",printID, false, 3, false)
writenewick_level1(net::HybridNetwork,outgroup::AbstractString) = writenewick_level1(net,false, true,true,outgroup,true, false, 3, false)
writenewick_level1(net::HybridNetwork,di::Bool,outgroup::AbstractString) = writenewick_level1(net,di, true,true,outgroup,true, false, 3, false)

"""
    writenewick_level1(net::HybridNetwork)

Write the extended Newick parenthetical format of a
level-1 network object with many optional arguments (see below).
Makes a deep copy of net: does *not* modify `net`.

- di=true: write in format for Dendroscope (default false)
- namelabel=true: If `namelabel` is true, taxa are labelled by their names;
otherwise taxa are labelled by their numbers (unique identifiers).
- outgroup (string): name of outgroup to root the tree/network.
  if "none" is given, the root is placed wherever possible.
- printID=true, only print branch lengths for identifiable egdes
  according to the snaq estimation procedure (default false)
  (true inside of `snaq!`.)
- round: rounds branch lengths and heritabilities γ (default: true)
- digits: digits after the decimal place for rounding (defult: 3)
- string: if true (default), returns a string,
  otherwise returns an IOBuffer object.
- multall: (default false). set to true when there are multiple
  alleles per population.

The topology may be written using a root different than net.rooti,
if net.rooti is incompatible with one of more hybrid node.
Missing hybrid names are written as "#Hi" where "i" is the hybrid node number if possible.
""" #"
writenewick_level1(net::HybridNetwork; di=false::Bool, string=true::Bool, namelabel=true::Bool,
    outgroup="none"::AbstractString, printID=false::Bool, round=false::Bool, digits=3::Integer,
    multall=false::Bool) =
writenewick_level1(net, di, string, namelabel, outgroup, printID, round, digits, multall)

# function to check if root is well-placed
# and look for a better place if not
# searches on net.node because net.rooti is the index in net.node
# if we search in net.edge, we then need to search in net.node
# this function is only used inside writenewick_level1
function updateRoot!(net::HybridNetwork, outgroup::AbstractString)
    checkroot = false
    if(outgroup == "none")
        @debug "no outgroup defined"
        checkroot = true
    else
        println("outgroup defined $(outgroup)")
        index = findfirst(n -> outgroup == n.name, net.node)
        index != nothing ||
            error("outgroup $(outgroup) not in net.names $(net.names)")
        node = net.node[index]
        node.leaf || error("outgroup $(outgroup) is not a leaf in net")
        length(net.node[index].edge) == 1 || error("strange leaf $(outgroup), node number $(net.node[index].number) with $(length(net.node[index].edge)) edges instead of 1")
        edge = net.node[index].edge[1]
        if edge.containroot
            DEBUGC && @debug "creating new node in the middle of the external edge $(edge.number) leading to outgroup $(node.number)"
            othernode = getOtherNode(edge,node)
            removeEdge!(othernode,edge)
            removeNode!(othernode,edge)
            max_edge = maximum([e.number for e in net.edge]);
            max_node = maximum([e.number for e in net.node]);
            newedge = Edge(max_edge+1) #fixit: maybe this edge not identifiable, need to add that check
            newnode = Node(max_node+1,false,false,[edge,newedge])
            if(cleaned(net) && !isTree(net) && !isempty(net.partition)) # fixit: this will crash if network estimated with snaq, and then manipulated
                part = whichpartition(net,edge)
                push!(net.partition[part].edges,newedge)
            end
            setNode!(edge,newnode)
            setNode!(newedge,newnode)
            setEdge!(othernode,newedge)
            setNode!(newedge,othernode)
            pushEdge!(net,newedge)
            pushNode!(net,newnode)
            t = edge.length
            if t == -1
                edge.length = -1
                newedge.length = -1
            else    
                setLength!(edge,t/2)
                setLength!(newedge,t/2)
            end
            net.rooti = length(net.node) #last node is root
       else
            @warn "external edge $(net.node[index].edge[1].number) leading to outgroup $(outgroup) cannot contain root, root placed wherever"
            checkroot = true
       end
    end
    if(checkroot && !isTree(net))
        checkRootPlace!(net)
    end
 end

# function to check if a node could be root
# by the containroot attribute of edges around it
function canBeRoot(n::Node)
    !n.hybrid || return false
    #!hasHybEdge(n) || return false #need to allow for some reason, check ipad notes
    !n.leaf || return false
    return any([e.containroot for e in n.edge])
end

# function to delete the extra node created in updateRoot
# this extra node is needed to be able to compare networks with the distance function
# but if left in the network, everything crashes (as everything assumes three edges per node)
# fromUpdateRoot=true if called after updateRoot (in which case leaf has to be a leaf), ow it is used in readTopUpdate
function undoRoot!(net::HybridNetwork, fromUpdateRoot::Bool=true)
    if length(getroot(net).edge) == 2
        root = getroot(net)
        leaf = getOtherNode(root.edge[1],root).leaf ? getOtherNode(root.edge[1],root) : getOtherNode(root.edge[2],root)
        (fromUpdateRoot && leaf.leaf) || error("root should have one neighbor leaf which has to be the outgroup defined")
        deleteIntLeafWhile!(net,root,leaf);
    end
end

# .out file from snaq written by optTopRuns
"""
    readsnaqnetwork(output file)

Read the estimated network from a `.out` file generated by [`snaq!`](@ref).
The network score is read also, and stored in the network's field `.loglik`.

Warning: despite the name "loglik", this score is only proportional
to the network's pseudo-deviance: the lower, the better.
Do NOT use this score to calculate an AIC or BIC (etc.) value.
"""
function readsnaqnetwork(file::AbstractString)
    open(file) do s
        line = readline(s)
        line[1] == '(' ||
          error("output file $(file) does not contain a tree in the first line, instead it has $(line); or we had trouble finding ploglik.")
        # println("Estimated network from file $(file): $(line)")
        net = readnewick(line)
        # readTopologyUpdate is inadequate: would replace missing branch lengths, which are unidentifiable, by 1.0 values
        try
            vec = split(line,"-Ploglik = ")
            loglik!(net, parse(Float64,vec[2]))
        catch e
            @warn "could not find the network score; the error was:"
            rethrow(e)
        end
        return net
    end
end

# function to change negative branch lengths to 1.0 for starting topology
# and to change big branch lengths to 10.0
# also uses setLength for all edges
function cleanBL!(net::HybridNetwork)
    ##println("missing branch lengths will be set to 1.0")
    for e in net.edge
        if(e.length < 0)
            setLength!(e,1.0)
        elseif(e.length > 10.0)
            setLength!(e,10.0)
        else
            setLength!(e,e.length)
        end
    end
end
