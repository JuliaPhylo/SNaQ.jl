using SNaQ, Random, StatsBase


function generate_tree(n::Int, seed::Int=42)
    n > 0 || error("n must be >0 (n=$(n))")
    rng = Random.seed!(seed)

    newick = "(" * join(["t$(j)" for j=1:n], ",") * ");"
    tre0 = readnewick(newick)
    for E in tre0.edge
        E.length = 0.0
    end
    tre = simulatecoalescent(tre0, 1, 1)[1]
    for node in tre.node
        if node.name == ""
            node.name = "i$(abs(rand(rng, Int)) % 10000)"
        end
    end
    return tre
end


function generate_net(n::Int, h::Int, seed::Int=42)
    (n > 0 && h > 0) || error("n and h must be >0 (n=$(n), h=$(h))")
    rng = Random.seed!(seed)

    net = generate_tree(n, seed)
    for j = 1:h
        SNaQ.addhybrid!(net, SNaQ.sampleaddhybridparameters(net, rng)...)
    end
    for node in net.node
        if node.name == ""
            node.name = "i$(abs(rand(rng, Int)) % 10000)"
        end
    end
    for E in net.edge
        if E.length < 0
            E.length = rand()
        end
    end
    for H in net.hybrid
        if getparentedge(H).gamma < 0 || getparentedgeminor(H).gamma < 0
            γ = rand()
            γ = min(γ, 1.0 - γ)
            getparentedge(H).gamma = 1.0 - γ
            getparentedgeminor(H).gamma = γ
        end
    end
    return net
end

# function opt_benchmark(ntaxa::Int=15, nhyb::Int=3, seed::Int=42)
#     Random.seed!(seed)
#     net = generate_net(ntaxa, nhyb, seed)
#     gts = simulatecoalescent(net, 1000, 1);
#     q, t = countquartetsintrees(gts, showprogressbar=false);
#     blocks, _ = findquartetequations(net);

#     return @benchmark optimize!(net, blocks, q)
# end

# function opt_profile(ntaxa::Int=15, nhyb::Int=3, seed::Int=42; nrep::Int=1)
#     Random.seed!(seed)
#     net = generate_net(ntaxa, nhyb, seed)
#     gts = simulatecoalescent(net, 1000, 1);
#     q, t = countquartetsintrees(gts, showprogressbar=false);
#     blocks, _ = findquartetequations(net);

#     start_time = time()
#     @profview for _ = 1:nrep optimize!(net, blocks, q) end
#     @info "Profiling took $(round(time() - start_time, digits=4)) seconds."
# end


# function eqn_profile(ntaxa::Int=15, nhyb::Int=3, seed::Int=42; nrep::Int=1)
#     Random.seed!(42)
#     net = generate_net(ntaxa, nhyb, seed)

#     start_time = time()
#     @profview for _ = 1:nrep findquartetequations(net) end
#     @info "Profiling took $(round(time() - start_time, digits=4)) seconds."
# end

