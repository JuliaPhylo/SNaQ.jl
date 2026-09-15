# Restricting network search spaces

## Overview

SNaQ was originally developed to infer binary, semi-directed, *level-1* phylogenetic
networks, but this level-1 restriction is no longer present. Now, by default, SNaQ
infers networks from the space of all binary, semi-directed phylogenetic networks,
and provides the functionality for the user to provide custom restrictions to the
network search space. The flexibility to explore different network spaces is controlled
via the `restrictions` keyword argument in [`snaq!`](@ref).

## What are restrictions?

Restrictions are functions that take a `HybridNetwork` as input and return `true` if the
network satisfies the restriction, and `false` otherwise. These functions are used to filter
which networks are considered during the search. Only networks that return `true` with the
provided restriction function will be considered as candidate networks during optimization.
Well-written topological restrictions are typically orders of magnitude faster than
optimization or likelihood computation functions, so using a priori knowledge to restrict
the search space can lead to improvements in inference runtime by ruling out candidate
networks that are known a priori to be poor fits.

## Built-in restriction functions

SNaQ provides several pre-built functions for common restriction scenarios.
These can be used individually or combined using the [`restrictionset`](@ref) function.

- [`restrictgalledtree`](@ref): Restrict to level-1 networks (galled trees).
- [`restrictgallednetwork`](@ref): Restrict to galled networks (allow higher levels).
- [`restrictmaximumlevel(level)`](@ref): Restrict the maximum level (number of hybrids in the most complex biconnected component).
- [`restrictrootedtreechild`](@ref): Restrict to rooted tree-child networks.
- [`restrictweaklytreechild`](@ref): Restrict to weakly tree-child networks.
- [`restrictstronglytreechild`](@ref): Restrict to strongly tree-child networks.
- [`norestrictions`](@ref): Allow any binary, semi-directed network.
- [`tcgidentifiable`](@ref): Restrict to tree-child galled networks with identifiable parameters.
- [`defaultrestrictions`](@ref): Equivalent to `norestrictions` — no restrictions are applied.

### Combining restrictions with `restrictionset`

The [`restrictionset`](@ref) function allows you to combine multiple restrictions into a single function.
It takes optional keyword arguments to specify which restrictions to apply:

```julia
snaq!(net0, d; restrictions=restrictionset(
    max_level=2,
    weakly_tree_child=true
))
```

Available keyword arguments for `restrictionset`:

- **`max_level::Real`** (default `Inf`): Maximum network level
- **`galled_tree::Bool`** (default `false`): Restrict to galled trees (level-1)
- **`galled_network::Bool`** (default `false`): Restrict to galled networks
- **`rooted_tree_child::Bool`** (default `false`): Restrict to rooted tree-child networks
- **`weakly_tree_child::Bool`** (default `false`): Restrict to weakly tree-child networks
- **`strongly_tree_child::Bool`** (default `false`): Restrict to strongly tree-child networks

If both `galled_tree` and `galled_network` are specified, `galled_tree` takes precedence.

## Defining custom restriction functions

If the built-in restrictions don't meet your needs, you can define custom restriction
functions. A restriction function must:

1. Take a single positional argument: a `HybridNetwork` object
2. Return `true` if the network satisfies your restriction, or `false` otherwise

### Example: known siblings or multiple samples

Suppose you have two samples for each taxa in your dataset, so you know that
these samples should always be sister to one another. If the samples are named
A1, A2, B1, B2, and so on, we could enforce this topological constraint
with the following restriction function:

```julia
using PhyloNetworks

function sistersiblings(network::HybridNetwork)
  for leaf in network.leaf
    # Does this leaf have an immediate sister?
    parentnode = getparent(leaf)
    siblings = getchildren(parentnode)

    # Here, one of the entries in the `siblings` variable will be
    # the leaf we're currently examining (`leaf`), while the other
    # is its sibling node. We use the ternary operator below to
    # select the sibling of the node `leaf`.
    sibling = siblings[1] == leaf ? siblings[2] : siblings[1]

    # If this sibling is an internal node rather than a leaf then
    # our restriction is not met
    if !sibling.leaf return false end

    # If the start of the sibling's name is not the same as the
    # start of `leaf`'s name then our restriction is not met
    prefix = sibling.name[1:(length(sibling.name)-1)]
    if !startswith(leaf.name, prefix) return false end
  end

  # If we got here, then none of the leaves in the network 
  # violated the restriction, so we can return true
  return true
end
snaq!(net0, d; restrictions=sistersiblings)
```

### Example: the set of hybrid taxa is known

Suppose you know exactly which taxa are hybrids. To implement this restriction,
we do **not** want to force every proposed network to have **all** of these
taxa placed as hybrids. It is highly unlikely that the search would find such
a network on its very first step, so we would get stuck in a loop where 0
proposed networks meet the restriction. Additionally, if the network has multiple
reticulations, we typically want to infer networks with progressively increasing
numbers of reticulations.

Instead, to implement this restriction we will only reject proposed networks
that have any hybrid nodes with descendants that are **not** a part of our known
set of hybrid taxa. We would do this as follows:

```julia
using PhyloNetworks
global known_hybrids = ["A", "B", "C", "D"];

function onlyknownhybrids(network::HybridNetwork)
  for H in network.hybrid
    # We can utilize an internal function of SNaQ to make
    # implementing this restriction a little easier:
    descendants = SNaQ.getleavesbelow([H])
    
    # If any of these descendants are not know hybrids,
    # our restriction is not met
    if any(d -> !(d.name in known_hybrids), descendants)
      return false
    end
  end

  # If we got here, then our restriction was not violated
  return true
end
snaq!(net0, d; restrictions=onlyknownhybrids)
```

Note that this approach also has the benefit of allowing all trees when `hmax=0`.

## Custom restrictions with parallel computing

When utilizing custom restrictions with parallel computing, the restriction
function needs to be defined on every worker process. This can be achieved by
defining the function within an `@everywhere` block:

```julia
using Distributed

@everywhere begin
  function customrestriction(network::HybridNetwork)::Bool
    # restriction code...
  end
end
```

If this step is not performed, anomalous errors from worker processes similar
to the following may be thrown:

```
LoadError: On worker 7:
UndefVarError: `#510#511` not defined in `SNaQ`
```
