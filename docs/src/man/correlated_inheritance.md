# Correlated inheritance

By default, SNaQ assumes independent inheritance: lineages coalescing into a hybrid node
independently choose which parental lineage to follow. This is controlled by the parameter `ρ` (rho), the
inheritance correlation parameter:

- `ρ = 0.0` (default): independent inheritance
- `ρ = 1.0`: completely dependent inheritance
- `0 < ρ < 1`: partially dependent inheritance

!!! info "Note:"
	The character ρ can be typed in the Julia REPL by typing `\rho` and then pressing the TAB key.

The correlated inheritance model is described in
[Fogg et al. 2023](https://doi.org/10.1093/sysbio/syad030).

## Network search with correlated inheritance

To search for a network under a correlated inheritance model, pass `ρ` to [`snaq!`](@ref):

```julia
net1 = snaq!(net0, dataCF, hmax=1, ρ=0.5, filename="net1_corr")
```

## Optimizing parameters with correlated inheritance

To optimize branch lengths and inheritance probabilities on a fixed topology under correlated
inheritance, pass `ρ` as the third positional argument to [`optimize!`](@ref):

```julia
optnet = deepcopy(truenet)
optimize!(optnet, raxmlCF, 0.5)  # composite log-likelihood: the higher, the better
loglik(optnet)                   # to access again later
```

## Computing the composite log-likelihood with correlated inheritance

[`computeloss`](@ref) also accepts `ρ` as its third positional argument:

```julia
computeloss(truenet, raxmlCF, 0.5)  # composite log-likelihood: the higher, the better
loglik(truenet)                     # to access again later
```

Note that this approach does not optimize the branch lengths or inheritance proportions of the provided
network, it only computes the composite log-likelihood (also called the loss) of the network with the
provided parameters.

## Conversion utilities

SNaQ internally uses a reparametrization `α = (1 - ρ) / ρ` (with `α = Inf` corresponding to `ρ = 0` and `α = 0` to `ρ = 1`).
Convenience functions are provided if you need to work with `α` directly:

```julia
α = rhotoalpha(0.5)  # ρ = 0.5 → α = 1.0
ρ = alphatorho(1.0)  # α = 1.0 → ρ = 0.5
```
