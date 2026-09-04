# LayeredMixedMatrices.jl

This Julia package is a proof-of-concept implementation of the combinatorial canonical
form (CCF) of a layered mixed matrix, presented in the forthcoming preprint
**Layered mixed matrices and reaction networks** by Arne Kuhrs, Máté L. Telek, and
Nicola Vassena.

The CCF is the combinatorially unique finest block-triangular form of a layered mixed
matrix `[T | Q]` under permutations and a change of basis of the numeric part, see
Murota, *Matrices and Matroids for Systems Analysis* (Springer, 2000). The
implementation follows the algorithm of Murota and Scharbrodt,
[*Computing the combinatorial canonical form of a layered mixed matrix*](https://doi.org/10.1080/10556789808805720),
Optimization Methods and Software, 1998.

For a reaction network with stoichiometric matrix `S` and symbolic reactivity matrix
`R`, the symbolic Jacobian is the product `SR`, and `det(SR) = c * det([R | C])` for a
Gale dual `C` of `S` and a nonzero constant `c`. The CCF of `[R | C]` then gives, in one
computation, the irreducible factorization of `det(SR)`, the buffering structures of
the network as the order ideals of the block poset, and the influence graph. The
package has an interface to [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and
[Oscar.jl](https://github.com/oscar-system/Oscar.jl).

## Installation
To install the package, run the following commands:
```julia
using Pkg
Pkg.add(url = "https://github.com/ArneKuhrs/LayeredMixedMatrices.jl")
```

## Examples of usage

You load the package in a Julia session by running the following command:

```julia
using LayeredMixedMatrices
```

You can either compute the CCF directly for a layered mixed matrix given by its
numeric and generic parts, or for a chemical reaction network given in Catalyst
format.

For the running example in the paper, a network with three species and six reactions,

```
X1 -> 2 X1                  X2 -> 2 X2                  X3 -> 2 X3
2 X1 + X2 -> X2             X2 -> X1 + 3 X2             X3 -> 0
```

we supply a basis `C` of `ker(S)` together with the support `Rpat` of the symbolic
reactivity matrix, where `Rpat[rho, m]` is `true` exactly when species `X_m` is a
reactant of reaction `rho`, and compute the CCF as follows.

```julia-repl
julia> C = [1 2 0; 0 1 0; 2 0 0; -1 0 0; 0 0 1; 0 0 1];

julia> Rpat = Bool[1 0 0; 1 1 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1];

julia> res = ccf(Rpat, C; columnLM = true)
CCFResult
  Columns:
    C0  = Int64[]
    C1 = [5, 6]
    C2 = [3, 4]
    C3 = [1, 2]
    C∞  = Int64[]
  Column-block order pairs = [(2, 3)]
  Column-block covers      = [(2, 3)]
  Q-row blocks:
    R0Q = Int64[]
    R1Q = [3]
    R2Q = [1]
    R3Q = [2]
    R∞Q = Int64[]
  T-row blocks:
    R0T = Int64[]
    R1T = [3]
    R2T = [2]
    R3T = [1]
    R∞T = Int64[]
  rowperm = [3, 6, 1, 5, 2, 4]
  colperm = [5, 6, 3, 4, 1, 2]
```

The three diagonal blocks are the subnetworks `({X1}, {1,2})`, `({X2}, {3,4})` and
`({X3}, {5,6})`, recovered from the network structure alone. Their determinants are
the irreducible factors of the symbolic Jacobian determinant.

The partial order on the blocks, whose order ideals are the buffering structures of
the network, is read off with the following commands.

```julia-repl
julia> block_poset(res)
1-element Vector{Tuple{Int64, Int64}}:
 (2, 3)

julia> column_block_order(res)
3×3 BitMatrix:
 1  0  0
 0  1  1
 0  0  1
```

We can also apply the same computation directly to a chemical reaction network,
through our interface to `Catalyst.jl` and `Oscar.jl` in the `examples` environment,
which additionally returns the symbolic determinants of the diagonal blocks.

```julia-repl
julia> using Pkg; Pkg.instantiate()

julia> using LayeredMixedMatrices, Oscar, Catalyst

julia> include("examples/oscar_interface.jl");

julia> rn = @reaction_network begin
           @species X1(t) X2(t) X3(t)
           k1, X1 --> 2*X1
           k2, 2*X1 + X2 --> X2
           k3, X2 --> 2*X2
           k4, X2 --> X1 + 3*X2
           k5, X3 --> 2*X3
           k6, X3 --> 0
       end;

julia> S = matrix(ZZ, netstoichmat(rn));

julia> Sminus = matrix(ZZ, substoichmat(rn));

julia> T, _, R = get_symbolic_matrix(transpose(Sminus));

julia> C = nullspace(S)[2];

julia> res = ccf(R, C; columnLM = true);

julia> for i in eachindex(res.Cblocks)
           B = vcat(matrix(T, res.P[res.RblocksQ[i], res.Cblocks[i]]),
                    transpose(R)[res.RblocksT[i], res.Cblocks[i]])
           println(i, ": ", det(B))
       end
1: -r_5_3 + r_6_3
2: r_3_2 + 2*r_4_2
3: -r_1_1 + 2*r_2_1

julia> factor(det(hcat(matrix(T, C), R)))
(r_3_2 + 2*r_4_2) * (r_1_1 - 2*r_2_1) * (r_5_3 - r_6_3)
```

To obtain the block-triangular matrix itself, we use the `ccf_form` command.

```julia-repl
julia> ccf_form(R, C; columnLM = true)
[2   r_1_1    1       0   0       0]
[1   r_2_1    0   r_2_2   0       0]
[0       0    2   r_3_2   0       0]
[0       0   -1   r_4_2   0       0]
[0       0    0       0   1   r_5_3]
[0       0    0       0   1   r_6_3]
```

## API

`ccf(A; maxiter)` implements the improved algorithm of Murota and Scharbrodt, whose
Steps A and B compute a large initial independent assignment before the
augmenting-path phase. `ccf_original(A; maxiter)` is Murota's original algorithm, starting from the empty assignment. Since the CCF is combinatorially unique the two must agree on every input, so it serves as a correctness check.

Note that the two modes take their arguments in opposite orders:

| call | first argument | second argument |
|---|---|---|
| `ccf(A, B; columnLM = true)`  | generic part | rational part |
| `ccf(A, B; columnLM = false)` | rational part | generic part |

A `CCFResult` carries the column blocks `C0, Cblocks, C∞`, the row blocks `R0Q/R0T`,
`RblocksQ/RblocksT`, `R∞Q/R∞T`, the pivoted matrix `P`, and the permutations
`rowperm`, `colperm`. It is read out with

| function | returns |
|---|---|
| `column_block_order(res)` | reflexive Boolean matrix of the block partial order |
| `column_block_order_pairs(res)` | the same order as pairs `(k,l)` |
| `block_poset(res)` | cover relations (Hasse edges) |
| `ccf_pattern_matrix(A, res)` | the block-triangular support pattern |
| `reverse_block_order(res)` | permutations listing blocks in reversed order |


## Reproducing the computations in the paper

From the repository root, with the `examples` environment instantiated as above:

```julia-repl
julia> include("examples/paper_computations.jl")
```

`examples/Manifest.toml` pins the exact package versions used. Two networks are treated: the running example above, and the *nominal cell model* of
Schliemann et al., with 47 species and 106 reactions, available as
[BIOMD0000000407](https://www.ebi.ac.uk/biomodels/BIOMD0000000407) in BioModels and
in [ODEbase](https://www.odebase.org/).

For the latter the CCF has 50 diagonal blocks and is computed in milliseconds. Of
these, 26 carry symbolic variables and hence give the 26 irreducible factors of the
symbolic Jacobian determinant:

| size of block | 2 | 2 | 2 | 3 | 5 | 5 | 7 | 13 | 26 |
|---|---|---|---|---|---|---|---|---|---|
| variables | 2 | 2 | 2 | 4 | 5 | 6 | 8 | 15 | 36 |
| monomials | 2 | 2 | 2 | 3 | 4 | 4 | 11 | 114 | 57341 |

(blocks of size 1 omitted). The expanded determinant therefore has 27,611,755,776
monomials. Computing all 50 block determinants takes a few seconds, whereas expanding
the full symbolic determinant directly does not terminate within 30 minutes.

## How to cite the package

If you find LayeredMixedMatrices.jl useful in your work, we kindly ask you to cite the
accompanying paper.

## License

MIT. See [LICENSE](LICENSE).