# LayeredMixedMatrices.jl

A Julia package computing the **combinatorial canonical form (CCF)** of a layered
mixed matrix, together with an interface for the symbolic Jacobians of chemical
reaction networks.

This is the companion code for the forthcoming preprint

> Arne Kuhrs, Máté L. Telek, Nicola Vassena.
> *Layered mixed matrices and reaction networks.*

## Background

A **layered mixed matrix** (LM-matrix) is a block matrix `[T | Q]` in which the
nonzero entries of the generic part `T` are algebraically independent over the field
to which the entries of the rational part `Q` belong. Its **combinatorial canonical
form** is the combinatorially unique finest block-triangular form attainable under
LM-equivalence, that is, under row and column permutations together with a change of
basis of the numeric part. The block structure is uniquely determined; the entries of
the numeric part are not. When the `Q`-part is empty, the CCF reduces to the
Dulmage–Mendelsohn decomposition.

The CCF was established by Murota, Iri and Nakamura; see Murota, *Matrices and
Matroids for Systems Analysis* (Springer, 2000) for a full account. It is computable
in polynomial time by matroid-theoretic methods, via a submodular flow model. Two of
its properties drive everything below: the determinant factors as the product of the
determinants of the diagonal blocks, and the block poset governs the zero pattern of
the inverse.

For a reaction network with stoichiometric matrix `S` and symbolic reactivity matrix
`R`, the symbolic Jacobian is `G = SR`. Choosing a Gale dual `C` of `S`, i.e. a basis
of `ker(S)`, the matrix `A = [R | C]` is an LM-matrix with `det(G) = c * det(A)` for
some nonzero constant `c`. The CCF of `A` then yields, in one computation:

* the irreducible factorization of `det(G)` as a polynomial in the reaction-rate
  derivatives,
* the buffering structures of the network, as the order ideals of the block poset,
* the influence graph, i.e. which reactions influence which species concentrations.

The implementation follows

> K. Murota and M. Scharbrodt. *Computing the combinatorial canonical form of a
> layered mixed matrix.* Optimization Methods and Software, 1998.
> [doi:10.1080/10556789808805720](https://doi.org/10.1080/10556789808805720)

The computation runs in two phases: a maximum independent assignment is found in an
auxiliary network, and that network is then decomposed into strongly connected
components.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/ArneKuhrs/LayeredMixedMatrices.jl")
```

## Example

We use the running example of the paper: three species, six reactions.

```
X1 -> 2 X1                  X2 -> 2 X2                  X3 -> 2 X3
2 X1 + X2 -> X2             X2 -> X1 + 3 X2             X3 -> 0
```

The matrix `C` is a basis of `ker(S)`, and `Rpat` is the support of the symbolic reactivity
matrix: `Rpat[rho, m]` is `true` exactly when species `X_m` is a reactant of
reaction `rho`.

```julia
using LayeredMixedMatrices

C = [ 1  2  0
      0  1  0
      2  0  0
     -1  0  0
      0  0  1
      0  0  1 ]

Rpat = Bool[ 1 0 0
             1 1 0
             0 1 0
             0 1 0
             0 0 1
             0 0 1 ]

res = ccf(Rpat, C; columnLM = true)
```

This prints the full decomposition. The essential part:

```
C1 = [5, 6]     C2 = [3, 4]     C3 = [1, 2]
Column-block covers = [(2, 3)]
```

Three diagonal blocks, each 2x2. Reading them against the network:

| block | reactions | species | determinant of the diagonal block |
|---|---|---|---|
| `C1` | 5, 6 | `X3` | `r_6_3 - r_5_3` |
| `C2` | 3, 4 | `X2` | `r_3_2 + 2 r_4_2` |
| `C3` | 1, 2 | `X1` | `2 r_2_1 - r_1_1` |

The algorithm recovers the influence graph of the reaction network from the network
structure alone, that is, from the stoichiometric matrix `S` and the support pattern
of `R`. No kinetics, no rate constants and no numerical values enter the computation.
The product of the three block determinants is, up to a nonzero constant, the full
symbolic Jacobian determinant.

`C0` and `C∞`, the column index sets of Murota's horizontal and vertical tails, are
both empty here; for a square matrix that is exactly the statement that it is
nonsingular.

The poset has the single cover relation `C2 -> C3`, with `C1` incomparable to both.
By the main theorem of the paper, the inverse is fully dense exactly on the three diagonal blocks and on
the `(2,3)` block, and identically zero elsewhere. The three-element poset has six
order ideals, and these are the six buffering structures of the network.

### With Catalyst and Oscar

In the `examples` environment the same computation runs from the network itself, and
returns the actual symbolic determinants rather than only the block structure.
The `include` path below is relative, so run this from the repository root.

```julia
using LayeredMixedMatrices, Oscar, Catalyst
include("examples/oscar_interface.jl")

rn = @reaction_network begin
    @species X1(t) X2(t) X3(t)
    k1, X1 --> 2*X1
    k2, 2*X1 + X2 --> X2
    k3, X2 --> 2*X2
    k4, X2 --> X1 + 3*X2
    k5, X3 --> 2*X3
    k6, X3 --> 0
end

S      = matrix(ZZ, netstoichmat(rn))
Sminus = matrix(ZZ, substoichmat(rn))
T, _, R = get_symbolic_matrix(transpose(Sminus))
C = nullspace(S)[2]

res = ccf(R, C; columnLM = true)

for i in eachindex(res.Cblocks)
    B = vcat(matrix(T, res.P[res.RblocksQ[i], res.Cblocks[i]]),
             transpose(R)[res.RblocksT[i], res.Cblocks[i]])
    println(i, ": ", det(B))
end
```

```
1: -r_5_3 + r_6_3
2: r_3_2 + 2*r_4_2
3: -r_1_1 + 2*r_2_1
```

which agrees with factoring the determinant directly:

```julia
A = hcat(matrix(T, C), R)
factor(det(A))
# (r_3_2 + 2*r_4_2) * (r_1_1 - 2*r_2_1) * (r_5_3 - r_6_3)
```

`ccf_form(R, C; columnLM = true)` returns the block-triangular matrix itself:

```
[2   r_1_1    1       0   0       0]
[1   r_2_1    0   r_2_2   0       0]
[0       0    2   r_3_2   0       0]
[0       0   -1   r_4_2   0       0]
[0       0    0       0   1   r_5_3]
[0       0    0       0   1   r_6_3]
```

### Argument order

The two modes take their arguments in opposite orders:

| call | first argument | second argument |
|---|---|---|
| `ccf(A, B; columnLM = true)`  | generic part | rational part |
| `ccf(A, B; columnLM = false)` | rational part | generic part |

## API

Both entry points compute the same combinatorial canonical form and return the same
`CCFResult`. They differ only in how the independent assignment is initialised, which
is where the speed difference comes from.

| function | algorithm |
|---|---|
| `ccf(A; maxiter)` | The improved algorithm of Murota and Scharbrodt (their Section 3). Step A greedily matches the rows of the generic part, so that the corresponding augmenting paths have length one and need no pivoting. Step B, in the *new algorithm* variant of their Section 3.2, runs over the rows of the numeric part and pivots on the first nonzero entry lying in a column that is not yet matched, which yields a large initial independent set. The augmenting-path phase then proceeds from there. This is the recommended entry point and the one used for the computations in the paper. |
| `ccf_original(A; maxiter)` | The original algorithm of Murota, as presented in Section 2 of Murota–Scharbrodt: the same augmenting-path phase, but started from the empty assignment. Kept as a reference implementation. |

Since the CCF is combinatorially unique, the two must return identical block
structures on every input, which makes `ccf_original` a correctness check on `ccf`.

A `CCFResult` carries the column blocks `C0, Cblocks, C∞`, the corresponding row
blocks `R0Q/R0T`, `RblocksQ/RblocksT`, `R∞Q/R∞T`, the pivoted matrix `P`, and the
permutations `rowperm`, `colperm`.

| function | returns |
|---|---|
| `column_block_order(res)` | reflexive Boolean matrix of the block partial order |
| `column_block_order_pairs(res)` | the same order as pairs `(k,l)` |
| `block_poset(res)` | cover relations (Hasse edges) |
| `ccf_pattern_matrix(A, res)` | the block-triangular support pattern |
| `reverse_block_order(res)` | permutations listing blocks in reversed order |
| `print_matrix(io, M)` | aligned printing |

`ccf_form(T, Q; columnLM)`, which returns the block-triangular matrix itself, is not
part of the package: it is defined in `examples/oscar_interface.jl` and requires Oscar.

## Reproducing the computations in the paper

`examples/` is a separate environment pinning the exact package versions used.

```
julia --project=examples
julia> include("examples/paper_computations.jl")
```

It additionally requires Oscar, Catalyst, Graphs and GraphPlot, all resolved by
`examples/Manifest.toml`. `examples/oscar_interface.jl` converts a Catalyst
reaction network into the symbolic matrices.

Two networks are treated: the running example above, and
[BIOMD0000000407](https://www.ebi.ac.uk/biomodels/BIOMD0000000407)
(Schliemann et al., TNF pro/anti-apoptosis; 47 species, 106 reactions).

For the latter the CCF has 50 diagonal blocks and is computed in milliseconds.
Of these, 26 carry symbolic variables and hence give the 26 irreducible factors of
the symbolic Jacobian determinant:

| size of block | 2 | 2 | 2 | 3 | 5 | 5 | 7 | 13 | 26 |
|---|---|---|---|---|---|---|---|---|---|
| variables | 2 | 2 | 2 | 4 | 5 | 6 | 8 | 15 | 36 |
| monomials | 2 | 2 | 2 | 3 | 4 | 4 | 11 | 114 | 57341 |

(blocks of size 1 omitted). The expanded determinant therefore has
27,611,755,776 monomials. Computing all 50 block determinants takes a few seconds,
whereas expanding the full symbolic determinant directly does not terminate within
30 minutes.

## Citing

If you use this code, please cite the paper above.

## License

MIT. See [LICENSE](LICENSE).