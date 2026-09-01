module LayeredMixedMatrices

using LinearAlgebra
using Graphs

export LMMatrix, CCFResult, ccf, ccf_original, ccf_pattern_matrix,
       ccf_matrix_any, print_matrix, column_block_order,
       column_block_order_pairs, block_poset, reverse_block_order

# ============================================================
# LM-matrix
# ============================================================

"""
LM-matrix A = [Q; T]

Q    : mQ × n numeric matrix over a field K
Tpat : mT × n Bool matrix for the zero/nonzero pattern of T
"""
struct LMMatrix{T}
    Q::Matrix{T}
    Tpat::BitMatrix
end

function LMMatrix(Q::AbstractMatrix{T}, Tpat::AbstractMatrix{Bool}) where {T}
    size(Q, 2) == size(Tpat, 2) ||
        throw(ArgumentError("Q and Tpat must have the same number of columns"))

    return LMMatrix{T}(Matrix(Q), BitMatrix(Tpat))
end

mQ(A::LMMatrix) = size(A.Q, 1)
mT(A::LMMatrix) = size(A.Tpat, 1)
ncols(A::LMMatrix) = size(A.Q, 2)

# ============================================================
# Vertex encoding
#
# Vertices of G~ are RT ∪ CQ ∪ C.
#
#   RT(i) = i                         for i = 1..mT
#   CQ(j) = mT + j                    for j = 1..n
#   C(j)  = mT + n + j                for j = 1..n
# ============================================================

struct VertexMap
    mT::Int
    n::Int
end

@inline rt(vmap::VertexMap, i::Int) = i
@inline cq(vmap::VertexMap, j::Int) = vmap.mT + j
@inline c(vmap::VertexMap, j::Int)  = vmap.mT + vmap.n + j

@inline is_rt(vmap::VertexMap, v::Int) =
    1 <= v <= vmap.mT

@inline is_cq(vmap::VertexMap, v::Int) =
    (vmap.mT + 1) <= v <= (vmap.mT + vmap.n)

@inline is_c(vmap::VertexMap, v::Int) =
    (vmap.mT + vmap.n + 1) <= v <= (vmap.mT + 2vmap.n)

@inline rt_index(vmap::VertexMap, v::Int) = v
@inline cq_index(vmap::VertexMap, v::Int) = v - vmap.mT
@inline c_index(vmap::VertexMap, v::Int)  = v - vmap.mT - vmap.n

@inline nv_total(vmap::VertexMap) = vmap.mT + 2vmap.n

# ============================================================
# State
# ============================================================

mutable struct CCFState{T}
    A::LMMatrix{T}
    vmap::VertexMap

    # base[h] ∈ {0,1,...,n}, h = 1..mQ
    base::Vector{Int}

    # Current transformed Q-block and row transformation matrix
    P::Matrix{T}
    S::Matrix{T}

    # Reoriented matching arcs M° : C -> (RT ∪ CQ)
    M_oriented::Set{Tuple{Int,Int}}
end

function CCFState(A::LMMatrix{T}) where {T}
    mq = mQ(A)

    return CCFState(
        A,
        VertexMap(mT(A), ncols(A)),
        zeros(Int, mq),
        copy(A.Q),
        Matrix{T}(I, mq, mq),
        Set{Tuple{Int,Int}}(),
    )
end

# ============================================================
# Step 2 data
# ============================================================

struct Step2Data
    I::BitSet
    J::BitSet
    Splus::BitSet
    SplusQ::BitSet
    Sminus::BitSet
end

# ============================================================
# Result
# ============================================================

struct CCFResult{T}
    P::Matrix{T}
    S::Matrix{T}
    base::Vector{Int}

    R0Q::Vector{Int}
    R0T::Vector{Int}
    RblocksQ::Vector{Vector{Int}}
    RblocksT::Vector{Vector{Int}}
    R∞Q::Vector{Int}
    R∞T::Vector{Int}

    C0::Vector{Int}
    Cblocks::Vector{Vector{Int}}
    C∞::Vector{Int}

    # Partial order among the proper column blocks C1,...,Cb.
    # block_order[k,l] == true iff Ck ⪯ Cl.  The relation is reflexive.
    block_order::BitMatrix

    # Cover relations (Hasse edges) of the same poset, as pairs (k,l).
    block_covers::Vector{Tuple{Int,Int}}

    rowperm::Vector{Int}
    colperm::Vector{Int}
end

# ============================================================
# Utility: current matched columns
# ============================================================

function matched_column_indices(state::CCFState)
    vmap = state.vmap
    cols = BitSet()

    for (u, _) in state.M_oriented
        if is_c(vmap, u)
            push!(cols, c_index(vmap, u))
        end
    end

    return cols
end

function matched_heads(state::CCFState)
    heads = BitSet()

    for (_, v) in state.M_oriented
        push!(heads, v)
    end

    return heads
end

# ============================================================
# Gaussian pivot
# ============================================================

"""
Perform the pivot update after assigning column j to row h.

Murota formula:

base[h] := j
w := 1/P[h,j]

P[k,l] := P[k,l] - w*P[k,j]*P[h,l]
S[k,l] := S[k,l] - w*P[k,j]*S[h,l]
P[k,j] := 0
"""
function pivot_assign_column!(state::CCFState{T}, h::Int, j::Int) where {T}
    mq, n = size(state.P)

    @assert 1 <= h <= mq
    @assert 1 <= j <= n
    @assert !iszero(state.P[h, j])

    state.base[h] = j

    w = inv(state.P[h, j])

    Ph = copy(@view state.P[h, :])
    Sh = copy(@view state.S[h, :])
    Pcolj = copy(@view state.P[:, j])

    @inbounds for k in 1:mq
        k == h && continue

        α = w * Pcolj[k]

        if !iszero(α)
            for l in 1:n
                if l != j
                    state.P[k, l] -= α * Ph[l]
                end
            end

            for l in 1:mq
                state.S[k, l] -= α * Sh[l]
            end

            state.P[k, j] = zero(T)
        end
    end

    return state
end

# ============================================================
# Upgraded initialization: Step A
# ============================================================

"""
Step A of Murota--Scharbrodt upgraded algorithm.

Compute a greedy maximal matching in the T-part graph RT -> C.

Each selected arc f_i -> x_j is inserted into M° as x_j -> f_i.
"""
function stepA_greedy_T_matching!(state::CCFState)
    A = state.A
    vmap = state.vmap

    matched_rows = falses(mT(A))
    matched_cols = falses(ncols(A))

    @inbounds for i in 1:mT(A)
        matched_rows[i] && continue

        for j in 1:ncols(A)
            if !matched_cols[j] && A.Tpat[i, j]
                push!(state.M_oriented, (c(vmap, j), rt(vmap, i)))
                matched_rows[i] = true
                matched_cols[j] = true
                break
            end
        end
    end

    return state
end

# ============================================================
# Upgraded initialization: Step B, new algorithm
# ============================================================

"""
Step B of the new algorithm.

For each Q-row h, choose the first currently unmatched column j with P[h,j] != 0.
Then add x_j -> x_jQ to M°, set base[h]=j, and pivot.

This avoids repeated augmenting paths of type CQ -> C.
"""
function stepB_greedy_Q_basis!(state::CCFState{T}) where {T}
    A = state.A
    vmap = state.vmap

    matched_cols = matched_column_indices(state)

    @inbounds for h in 1:mQ(A)
        state.base[h] != 0 && continue

        chosen = 0

        for j in 1:ncols(A)
            if !(j in matched_cols) && !iszero(state.P[h, j])
                chosen = j
                break
            end
        end

        if chosen != 0
            # Add the reversed EQ arc x_j -> x_jQ to M°
            push!(state.M_oriented, (c(vmap, chosen), cq(vmap, chosen)))
            push!(matched_cols, chosen)

            # Assign base and pivot
            pivot_assign_column!(state, h, chosen)
        end
    end

    return state
end

# ============================================================
# Step 2
# ============================================================

function step2_sets(state::CCFState{T}) where {T}
    A = state.A
    vmap = state.vmap
    n = ncols(A)

    matched_left = BitSet()
    matched_cols = BitSet()

    for (u, v) in state.M_oriented
        push!(matched_left, v)
        push!(matched_cols, c_index(vmap, u))
    end

    # I = {j | jQ ∈ ∂⁻M° ∩ CQ}
    I = BitSet()

    @inbounds for j in 1:n
        if cq(vmap, j) in matched_left
            push!(I, j)
        end
    end

    # J = {j ∈ C \ I | all free base rows have P[h,j] = 0}
    free_rows = findall(==(0), state.base)

    J = BitSet()

    @inbounds for j in 1:n
        if j in I
            continue
        end

        ok = true

        for h in free_rows
            if !iszero(state.P[h, j])
                ok = false
                break
            end
        end

        if ok
            push!(J, j)
        end
    end

    # S_T^+
    Splus = BitSet()
    SplusQ = BitSet()

    @inbounds for i in 1:mT(A)
        v = rt(vmap, i)

        if !(v in matched_left)
            push!(Splus, v)
        end
    end

    # S_Q^+
    @inbounds for j in 1:n
        if !(j in I) && !(j in J)
            v = cq(vmap, j)
            push!(Splus, v)
            push!(SplusQ, v)
        end
    end

    # S^-
    Sminus = BitSet()

    @inbounds for j in 1:n
        if !(j in matched_cols)
            push!(Sminus, c(vmap, j))
        end
    end

    return Step2Data(I, J, Splus, SplusQ, Sminus)
end

# ============================================================
# Graph construction
# ============================================================

"""
Build G~ = (V~, E~).

E~ = ET ∪ EQ ∪ E+ ∪ M°

Important correction:
E+ is restricted to j ∈ J, as in Murota's definition.
"""
function build_graph(state::CCFState{T}, data::Step2Data) where {T}
    A = state.A
    vmap = state.vmap

    G = SimpleDiGraph(nv_total(vmap))

    # ET: RT -> C
    @inbounds for i in 1:mT(A), j in 1:ncols(A)
        if A.Tpat[i, j]
            Graphs.add_edge!(G, rt(vmap, i), c(vmap, j))
        end
    end

    # EQ: CQ -> C
    @inbounds for j in 1:ncols(A)
        Graphs.add_edge!(G, cq(vmap, j), c(vmap, j))
    end

    # E+: CQ(i) -> CQ(j), only for j ∈ J
    @inbounds for h in 1:mQ(A)
        i = state.base[h]

        if i != 0
            ui = cq(vmap, i)

            for j in data.J
                if !iszero(state.P[h, j])
                    Graphs.add_edge!(G, ui, cq(vmap, j))
                end
            end
        end
    end

    # M°
    for (u, v) in state.M_oriented
        Graphs.add_edge!(G, u, v)
    end

    return G
end

# ============================================================
# Multi-source BFS shortest path
# ============================================================

function shortest_path_any(G::SimpleDiGraph, sources::BitSet, sinks::BitSet)
    isempty(sources) && return nothing
    isempty(sinks) && return nothing

    nV = nv(G)

    seen = falses(nV)
    pred = zeros(Int, nV)

    q = Vector{Int}(undef, nV)
    qtail = 0

    for s in sources
        qtail += 1
        q[qtail] = s
        seen[s] = true
    end

    head = 1
    target = 0

    while head <= qtail
        u = q[head]
        head += 1

        if u in sinks
            target = u
            break
        end

        for v in outneighbors(G, u)
            if !seen[v]
                seen[v] = true
                pred[v] = u
                qtail += 1
                q[qtail] = v
            end
        end
    end

    target == 0 && return nothing

    path = [target]
    cur = target

    while pred[cur] != 0
        cur = pred[cur]
        pushfirst!(path, cur)
    end

    return path
end

# ============================================================
# Step 3
# ============================================================

function update_matching_along_path!(state::CCFState, path::Vector{Int})
    vmap = state.vmap

    # Remove M° arcs on path
    for k in 1:length(path)-1
        e = (path[k], path[k+1])

        if e in state.M_oriented
            delete!(state.M_oriented, e)
        end
    end

    # Reverse ET and EQ arcs into M°
    for k in 1:length(path)-1
        u = path[k]
        v = path[k+1]

        if is_rt(vmap, u) && is_c(vmap, v)
            push!(state.M_oriented, (v, u))
        elseif is_cq(vmap, u) && is_c(vmap, v)
            push!(state.M_oriented, (v, u))
        end
    end

    return state
end

function step3!(state::CCFState{T}, data::Step2Data, path::Vector{Int}) where {T}
    vmap = state.vmap

    update_matching_along_path!(state, path)

    # If the initial vertex belongs to S_Q^+, pivot it first
    startv = first(path)

    if startv in data.SplusQ
        j = cq_index(vmap, startv)

        h = findfirst(
            h -> state.base[h] == 0 && !iszero(state.P[h, j]),
            eachindex(state.base),
        )

        h === nothing && error("No admissible pivot row for initial S_Q^+ vertex.")

        pivot_assign_column!(state, h, j)
    end

    # For all E+ arcs CQ(i) -> CQ(j) on the path, pivot base i -> j
    for k in 1:length(path)-1
        u = path[k]
        v = path[k+1]

        if is_cq(vmap, u) && is_cq(vmap, v)
            i = cq_index(vmap, u)
            j = cq_index(vmap, v)

            h = findfirst(==(i), state.base)

            h === nothing && error("No row h with base[h] = $i.")

            pivot_assign_column!(state, h, j)
        end
    end

    return state
end

# ============================================================
# Reachability
# ============================================================

function reachable_from(G::SimpleDiGraph, starts::BitSet)
    nV = nv(G)

    seen = falses(nV)
    q = Vector{Int}(undef, nV)
    qtail = 0

    for s in starts
        qtail += 1
        q[qtail] = s
        seen[s] = true
    end

    head = 1

    while head <= qtail
        u = q[head]
        head += 1

        for v in outneighbors(G, u)
            if !seen[v]
                seen[v] = true
                qtail += 1
                q[qtail] = v
            end
        end
    end

    return BitSet(findall(seen))
end

function reachable_to(G::SimpleDiGraph, targets::BitSet)
    nV = nv(G)

    seen = falses(nV)
    q = Vector{Int}(undef, nV)
    qtail = 0

    for t in targets
        qtail += 1
        q[qtail] = t
        seen[t] = true
    end

    head = 1

    while head <= qtail
        u = q[head]
        head += 1

        for v in inneighbors(G, u)
            if !seen[v]
                seen[v] = true
                qtail += 1
                q[qtail] = v
            end
        end
    end

    return BitSet(findall(seen))
end

# ============================================================
# Step 4 helpers
# ============================================================

function induced_subgraph_with_maps(G::SimpleDiGraph, keep::Vector{Int})
    old_to_new = Dict{Int,Int}()

    for (i, v) in enumerate(keep)
        old_to_new[v] = i
    end

    H = SimpleDiGraph(length(keep))
    keep_set = Set(keep)

    for u in keep
        u2 = old_to_new[u]

        for v in outneighbors(G, u)
            if v in keep_set
                Graphs.add_edge!(H, u2, old_to_new[v])
            end
        end
    end

    return H, old_to_new, keep
end

function scc_condensation_dag(H::SimpleDiGraph, comps::Vector{Vector{Int}})
    cid = zeros(Int, nv(H))

    for (k, comp) in enumerate(comps), v in comp
        cid[v] = k
    end

    dag = SimpleDiGraph(length(comps))

    for e in edges(H)
        cu = cid[src(e)]
        cv = cid[dst(e)]

        if cu != cv && !has_edge(dag, cu, cv)
            Graphs.add_edge!(dag, cu, cv)
        end
    end

    return dag
end

function topo_order_scc_dag(H::SimpleDiGraph, comps::Vector{Vector{Int}})
    return topological_sort(scc_condensation_dag(H, comps))
end

# Reachability order restricted to those SCCs that actually contain column
# vertices.  This is exactly the order used in Step 4 of Murota's CCF
# algorithm: V_lambda ⪯ V_mu iff there is a directed path V_lambda -> V_mu.
function column_block_poset(dag::SimpleDiGraph, block_cids::Vector{Int})
    b = length(block_cids)
    order = falses(b, b)

    for k in 1:b
        reached = reachable_from(dag, BitSet([block_cids[k]]))
        for l in 1:b
            order[k, l] = block_cids[l] in reached
        end
    end

    covers = Tuple{Int,Int}[]
    for k in 1:b, l in 1:b
        k == l && continue
        order[k, l] || continue

        iscover = true
        for h in 1:b
            (h == k || h == l) && continue
            if order[k, h] && order[h, l]
                iscover = false
                break
            end
        end

        iscover && push!(covers, (k, l))
    end

    return BitMatrix(order), covers
end

# ============================================================
# Step 4
# ============================================================

function step4(state::CCFState{T}) where {T}
    A = state.A
    vmap = state.vmap

    data = step2_sets(state)
    G = build_graph(state, data)

    V∞ = reachable_from(G, data.Splus)
    V0 = reachable_to(G, data.Sminus)

    C0 = sort([c_index(vmap, v) for v in V0 if is_c(vmap, v)])
    C∞ = sort([c_index(vmap, v) for v in V∞ if is_c(vmap, v)])

    keep = sort([v for v in vertices(G) if !(v in V0) && !(v in V∞)])

    H, _, new_to_old = induced_subgraph_with_maps(G, keep)

    comps_local = strongly_connected_components(H)
    dag = scc_condensation_dag(H, comps_local)
    topo = topological_sort(dag)

    Cblocks = Vector{Vector{Int}}()
    comp_cols = Vector{BitSet}()
    comp_trows = Vector{Vector{Int}}()
    block_cids = Int[]

    for cid in topo
        comp_old = [new_to_old[v] for v in comps_local[cid]]

        cols = sort([c_index(vmap, v) for v in comp_old if is_c(vmap, v)])

        isempty(cols) && continue

        trows = sort([rt_index(vmap, v) for v in comp_old if is_rt(vmap, v)])

        push!(Cblocks, cols)
        push!(comp_cols, BitSet(cols))
        push!(comp_trows, trows)
        push!(block_cids, cid)
    end

    # Partial order among C1,...,Cb, induced by reachability among the SCCs.
    # Paths are allowed to pass through SCCs containing no column vertices.
    block_order, block_covers = column_block_poset(dag, block_cids)

    R0T = sort([rt_index(vmap, v) for v in V0 if is_rt(vmap, v)])
    R∞T = sort([rt_index(vmap, v) for v in V∞ if is_rt(vmap, v)])

    C0set = BitSet(C0)
    C∞set = BitSet(C∞)

    R0Q = sort([h for h in 1:mQ(A) if state.base[h] in C0set])
    R∞Q = sort([h for h in 1:mQ(A) if state.base[h] == 0 || state.base[h] in C∞set])

    RblocksQ = Vector{Vector{Int}}()
    RblocksT = Vector{Vector{Int}}()

    for (colset, trows) in zip(comp_cols, comp_trows)
        qrows = sort([h for h in 1:mQ(A) if state.base[h] in colset])

        push!(RblocksQ, qrows)
        push!(RblocksT, trows)
    end

    rowperm = Int[]

    append!(rowperm, R0Q)
    append!(rowperm, [mQ(A) + i for i in R0T])

    for b in eachindex(RblocksQ)
        append!(rowperm, RblocksQ[b])
        append!(rowperm, [mQ(A) + i for i in RblocksT[b]])
    end

    append!(rowperm, R∞Q)
    append!(rowperm, [mQ(A) + i for i in R∞T])

    colperm = Int[]

    append!(colperm, C0)

    for cols in Cblocks
        append!(colperm, cols)
    end

    append!(colperm, C∞)

    return CCFResult(
        copy(state.P),
        copy(state.S),
        copy(state.base),
        R0Q,
        R0T,
        RblocksQ,
        RblocksT,
        R∞Q,
        R∞T,
        C0,
        Cblocks,
        C∞,
        block_order,
        block_covers,
        rowperm,
        colperm,
    )
end

# ============================================================
# Main algorithms
# ============================================================

"""
Original algorithm: starts from empty M°.
"""
function ccf_original(A::LMMatrix{T}; maxiter::Int = 100_000) where {T}
    state = CCFState(A)

    for _iter in 1:maxiter
        data = step2_sets(state)
        G = build_graph(state, data)

        path = shortest_path_any(G, data.Splus, data.Sminus)

        if path === nothing
            return step4(state)
        end

        step3!(state, data, path)
    end

    error("Maximum iterations exceeded.")
end

"""
Upgraded algorithm.

Initialization:
- Step A: greedy matching in T-part.
- Step B: greedy initial Q-basis.

Then original augmenting-path phase.
"""
function ccf(A::LMMatrix{T}; maxiter::Int = 100_000) where {T}
    state = CCFState(A)

    # Upgraded initialization
    stepA_greedy_T_matching!(state)
    stepB_greedy_Q_basis!(state)

    # Original augmenting phase
    for _iter in 1:maxiter
        data = step2_sets(state)
        G = build_graph(state, data)

        path = shortest_path_any(G, data.Splus, data.Sminus)

        if path === nothing
            return step4(state)
        end

        step3!(state, data, path)
    end

    error("Maximum iterations exceeded.")
end


function ccf(Q, T; columnLM::Bool = true, kwargs...)
    if columnLM
        Qmat = Rational{BigInt}.(Array(transpose(T)))
        Tpat = BitMatrix(Array(transpose(Q)) .!= 0)
        return ccf(LMMatrix(Qmat, Tpat); kwargs...)
    else
        Qmat = Rational{BigInt}.(Array(Q))
        Tpat = BitMatrix(Array(T) .!= 0)
        return ccf(LMMatrix(Qmat, Tpat); kwargs...)
    end
end


# ============================================================
# Column-block partial order helpers
# ============================================================

"""
    column_block_order(res::CCFResult)

Return the reflexive Boolean matrix of the partial order on the proper
column blocks `C1,...,Cb`.  Entry `(k,l)` is `true` exactly when
`Ck ⪯ Cl`, i.e. when the SCC defining `Ck` can reach the SCC defining `Cl`
in the final graph of Step 4.
"""
column_block_order(res::CCFResult) = copy(res.block_order)

"""
    column_block_order_pairs(res::CCFResult; reflexive=false)

Return the partial order as pairs `(k,l)` of block indices.  By default the
reflexive pairs `(k,k)` are omitted.
"""
function column_block_order_pairs(res::CCFResult; reflexive::Bool = false)
    b = length(res.Cblocks)
    return [(k,l) for k in 1:b for l in 1:b
            if res.block_order[k,l] && (reflexive || k != l)]
end

"""
    block_poset(res::CCFResult)

Return the cover relations `(k,l)` (the Hasse-diagram edges) of the partial
order on `C1,...,Cb`.
"""
block_poset(res::CCFResult) = copy(res.block_covers)

"""
    reverse_block_order(res::CCFResult)

Return row and column permutations listing the blocks in the reversed order
`∞, b, …, 1, 0`.  Transposing a row-form CCF gives a *lower* block-triangular matrix,
so the block order has to be reversed to recover the upper block-triangular column form.
"""
function reverse_block_order(res::CCFResult)
    mq = size(res.P, 1)

    rowblocks = Vector{Vector{Int}}()
    push!(rowblocks, vcat(res.R0Q, mq .+ res.R0T))
    for k in eachindex(res.RblocksQ)
        push!(rowblocks, vcat(res.RblocksQ[k], mq .+ res.RblocksT[k]))
    end
    push!(rowblocks, vcat(res.R∞Q, mq .+ res.R∞T))

    colblocks = Vector{Vector{Int}}()
    push!(colblocks, res.C0)
    for Ck in res.Cblocks
        push!(colblocks, Ck)
    end
    push!(colblocks, res.C∞)

    return vcat(reverse(rowblocks)...), vcat(reverse(colblocks)...)
end

# ============================================================
# Display helpers
# ============================================================

function Base.show(io::IO, res::CCFResult)
    println(io, "CCFResult")

    println(io, "  Columns:")
    println(io, "    C0  = ", res.C0)

    for (k, Ck) in enumerate(res.Cblocks)
        println(io, "    C$k = ", Ck)
    end

    println(io, "    C∞  = ", res.C∞)
    println(io, "  Column-block order pairs = ", column_block_order_pairs(res))
    println(io, "  Column-block covers      = ", res.block_covers)

    println(io, "  Q-row blocks:")
    println(io, "    R0Q = ", res.R0Q)

    for (k, Rk) in enumerate(res.RblocksQ)
        println(io, "    R$(k)Q = ", Rk)
    end

    println(io, "    R∞Q = ", res.R∞Q)

    println(io, "  T-row blocks:")
    println(io, "    R0T = ", res.R0T)

    for (k, Rk) in enumerate(res.RblocksT)
        println(io, "    R$(k)T = ", Rk)
    end

    println(io, "    R∞T = ", res.R∞T)

    println(io, "  rowperm = ", res.rowperm)
    print(io, "  colperm = ", res.colperm)
end

function ccf_pattern_matrix(A::LMMatrix, res::CCFResult)
    mq = mQ(A)
    mt = mT(A)
    n = ncols(A)

    M = Matrix{String}(undef, mq + mt, n)

    for i in 1:mq, j in 1:n
        M[i, j] = iszero(res.P[i, j]) ? "0" : string(res.P[i, j])
    end

    for i in 1:mt, j in 1:n
        M[mq + i, j] = A.Tpat[i, j] ? "*" : "0"
    end

    return M[res.rowperm, res.colperm]
end

function ccf_matrix_any(res::CCFResult, Tblock::AbstractMatrix)
    Pany = Matrix{Any}(res.P)
    Tany = Matrix{Any}(Tblock)
    PT = vcat(Pany, Tany)

    return PT[res.rowperm, res.colperm]
end

function pretty_entry(x)
    if x isa AbstractString
        return x
    elseif x isa AbstractFloat
        if iszero(x)
            return "0"
        elseif isinteger(x)
            return string(Int(round(x)))
        else
            return string(x)
        end
    elseif x isa Integer
        return string(x)
    elseif x isa Rational
        return string(x)
    else
        return string(x)
    end
end

function print_matrix(io::IO, M::AbstractMatrix)
    m, n = size(M)

    S = Matrix{String}(undef, m, n)

    for i in 1:m, j in 1:n
        S[i, j] = pretty_entry(M[i, j])
    end

    widths = [
        maximum(length(S[i, j]) for i in 1:m)
        for j in 1:n
    ]

    for i in 1:m
        print(io, "[")

        for j in 1:n
            print(io, lpad(S[i, j], widths[j]))

            if j < n
                print(io, "   ")
            end
        end

        print(io, "]")

        if i < m
            println(io)
        end
    end
end

print_matrix(M::AbstractMatrix) = print_matrix(stdout, M)

end


