function get_symbolic_matrix(A)
    m, n = size(A)

    names = ["r_$(i)_$(j)" for i in 1:m for j in 1:n if !iszero(A[i, j])]
    T, vars = polynomial_ring(QQ, names)

    R = zero_matrix(T, m, n)

    k = 1
    for i in 1:m, j in 1:n
        if !iszero(A[i, j])
            R[i, j] = vars[k]
            k += 1
        end
    end

    return T, vars, R
end

function ccf_form(T::MatElem,Q::MatElem;columnLM::Bool = true)
    if  columnLM
        ccf_res = ccf(T,Q,columnLM=true)
        reversed_rowperm,reversed_colperm = reverse_block_order(ccf_res)
        return transpose((vcat(matrix(base_ring(T),ccf_res.P),transpose(T)))[reversed_rowperm,reversed_colperm])
    else
        #Here we change the role of Q and T to make the code consistent
        ccf_res = ccf(T,Q,columnLM=false)
        return (vcat(matrix(base_ring(Q),ccf_res.P),Q))[ccf_res.rowperm,ccf_res.colperm]
    end
end

function plot_column_block_order(res)
    pairs = column_block_order_pairs(res)
    b = length(res.Cblocks)

    G = SimpleDiGraph(b)

    for (k, l) in pairs
        add_edge!(G, k, l)
    end

    labels = ["C$i" for i in 1:b]

    gplot(
        G,
        nodelabel = labels,
        nodesize = 0.15,
        arrowlengthfrac = 0.08
    )
end

function plot_column_block_poset(res)
    pairs = block_poset(res)
    b = length(res.Cblocks)

    G = Graphs.SimpleDiGraph(b)

    for (k, l) in pairs
        Graphs.add_edge!(G, k, l)
    end

    labels = ["C$i" for i in 1:b]

    GraphPlot.gplot(
        G,
        nodelabel = labels,
        nodesize = 0.15,
        arrowlengthfrac = 0.08
    )
end