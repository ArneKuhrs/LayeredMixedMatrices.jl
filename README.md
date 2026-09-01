# LayeredMixedMatrices.jl

A Julia package computing the **combinatorial canonical form (CCF)** of a layered
mixed matrix, together with an interface for the symbolic Jacobians of chemical
reaction networks.

This is the companion code for the forthcoming preprint 

> A. Kuhrs, M. L. Telek, N. Vassena.
> *Layered mixed matrices and reaction networks.*

## Background

A **layered mixed matrix** (LM-matrix) is a matrix split into a rational part `Q`
and a generic part `T` whose nonzero entries are algebraically independent.
Murota's combinatorial canonical form is the finest block-triangular
decomposition attainable under LM-equivalence; it is canonical and computable in
polynomial time. The determinant factors as the product of the determinants of
the diagonal blocks, and the block poset governs the zero pattern of the inverse.

The implementation follows the augmenting-path algorithm of
Murota and Scharbrodt.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ArneKuhrs/LayeredMixedMatrices.jl")
```

## Quick start

No external dependencies are needed for the core routines:

```julia
using LayeredMixedMatrices

# rational part: a basis of ker(S)
C = [ 1  2  0
      0  1  0
      2  0  0
     -1  0  0
      0  0  1
      0  0  1 ]

# support of the generic part (true = algebraically independent entry)
Rpat = Bool[ 1 0 0
             1 1 0
             0 1 0
             0 1 0
             0 0 1
             0 0 1 ]

res = ccf(Rpat, C; columnLM = true)
```

`res.Cblocks` gives the three column blocks, and `block_poset(res)` their Hasse
diagram.

### Argument order

The two modes take their arguments in opposite orders:

| call | first argument | second argument |
|---|---|---|
| `ccf(A, B; columnLM = true)`  | generic part | rational part |
| `ccf(A, B; columnLM = false)` | rational part | generic part |

## API

`ccf(A; maxiter)` — CCF with greedy initialisation (Steps A and B).
`ccf_original(A; maxiter)` — the published algorithm, starting from an empty
matching. Both return the same `CCFResult`; the first is faster.

A `CCFResult` carries the column blocks `C0, Cblocks, C∞`, the corresponding row
blocks `R0Q/R0T`, `RblocksQ/RblocksT`, `R∞Q/R∞T`, the pivoted matrix `P`, and the
permutations `rowperm`, `colperm`.

| function | returns |
|---|---|
| `column_block_order(res)` | reflexive Boolean matrix of the block partial order |
| `column_block_order_pairs(res)` | the same order as pairs `(k,l)` |
| `block_poset(res)` | cover relations (Hasse edges) |
| `reverse_block_order(res)` | permutations listing blocks in reversed order |
| `ccf_pattern_matrix(A, res)` | the block-triangular support pattern |
| `print_matrix(io, M)` | aligned printing |

## Reproducing the computations in the paper

`examples/` is a separate environment pinning the exact package versions used.

```
julia --project=examples
julia> include("examples/paper_computations.jl")
```

It additionally requires Oscar, Catalyst, Graphs and GraphPlot, all resolved by
`examples/Manifest.toml`. `examples/oscar_interface.jl` converts a Catalyst
reaction network into the symbolic matrices.

Two networks are treated: the running example of the paper, and
[BIOMD0000000407](https://www.ebi.ac.uk/biomodels/BIOMD0000000407)
(Schliemann et al., TNF pro/anti-apoptosis; 47 species, 106 reactions).

For the latter the CCF has 50 diagonal blocks. Computing all 50 block
determinants takes a few seconds and yields 27,611,755,776 monomials in total,
whereas expanding the full symbolic determinant directly does not terminate
within 30 minutes.

## Citing

If you use this code, please cite the paper above.

## License

MIT. See [LICENSE](LICENSE).
