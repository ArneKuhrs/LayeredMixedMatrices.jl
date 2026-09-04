using LayeredMixedMatrices
using Oscar
using Catalyst
using Graphs, GraphPlot
include(joinpath(@__DIR__, "oscar_interface.jl"))

rn = @reaction_network begin
    @species X1(t) X2(t) X3(t)
        k1, X1 --> 2*X1
        k2, 2*X1  + X2 -->   X2
        k3, X2 --> 2*X2
        k4, X2 --> X1 +3*X2 
        k5, X3 --> 2*X3
        k6, X3 --> 0
end
#Set up the data
S = matrix(ZZ, netstoichmat(rn)) #stoichiometric matrix
Sminus = matrix(ZZ, substoichmat(rn)) #reactant matrix
n,r = size(S)
T, _, R = get_symbolic_matrix(transpose(Sminus)) #T is the polynomial ring where the symbols in Rt live in
C = nullspace(S)[2]
A = hcat(matrix(T,C),R)
factor(det(A))

#Row LM-matrices
MT = transpose(R)
MQ = transpose(C)
ccf_res_row = ccf(MQ,MT,columnLM = false)
L_ccf = ccf_form(MQ,MT,columnLM = false)

#Column CCF
ccf_res_col = ccf(R,C,columnLM = true)

for i in eachindex(ccf_res_col.Cblocks)
    B = vcat(matrix(T, ccf_res_col.P[ccf_res_col.RblocksQ[i], ccf_res_col.Cblocks[i]]),
             transpose(R)[ccf_res_col.RblocksT[i], ccf_res_col.Cblocks[i]])
    println(i, ": ", det(B))
end

reverse_block_order(ccf_res_col)
A_ccf = ccf_form(R,C,columnLM = true)
column_block_order(ccf_res_col)
column_block_order_pairs(ccf_res_col)
block_poset(ccf_res_col)

#BIG example BIOMD0000000407
#Schliemann2011_TNF_ProAntiApoptosis
rn = @reaction_network begin 
k1, TNFR --> TNFR_E
k2, 0 --> TNFR
k3, TNFR_E --> 0
(k4,l4), 0 <--> RIP
(k5,l5), 0 <--> TRADD
(k6,l6), 0 <--> TRAF2
(k7,l7), 0 <--> FADD
k8, TNF_TNFR_E --> 0
k9, TNF_TNFR_TRADD --> 0
k10, TNFRC1 --> 0
k11, TNFRC2 --> 0
k12, TNFRC2_FLIP --> 0
k13, TNFRC2_FLIP_FLIP --> 0
k14, TNFRC2_pCasp8 --> 0
k15, TNFRC2_pCasp8_pCasp8 --> 0
k16, TNFRC2_FLIP_pCasp8 --> 0
k17, TNFRC2_FLIP_pCasp8_RIP_TRAF2 --> 0
(k18,l18), TNFR_E + TNF_E <--> TNF_TNFR_E
k19, TNF_TNFR_E + TRADD --> TNF_TNFR_TRADD
k20, RIP + TNF_TNFR_TRADD + TRAF2 --> TNFRC1
k21, TNFRC1 --> TNFRCint1
k22, TNFRCint1 --> RIP + TNFRCint2 + TRAF2
k23, 2*FADD + TNFRCint2 --> TNFRCint3
k24, TNFRCint3 --> TNFRC2
k25, FLIP + TNFRC2 --> TNFRC2_FLIP
k26, FLIP + TNFRC2_FLIP --> TNFRC2_FLIP_FLIP
k27, TNFRC2 + pCasp8 --> TNFRC2_pCasp8
k28, TNFRC2_pCasp8 + pCasp8 --> TNFRC2_pCasp8_pCasp8
k29, TNFRC2_pCasp8_pCasp8 --> Casp8 + TNFRC2
k30, FLIP + TNFRC2_pCasp8 --> TNFRC2_FLIP_pCasp8
k31, TNFRC2_FLIP + pCasp8 --> TNFRC2_FLIP_pCasp8
k32, TNFRC2_FLIP_pCasp8 --> Casp8 + TNFRC2
k33, RIP + TNFRC2_FLIP_pCasp8 + TRAF2 --> TNFRC2_FLIP_pCasp8_RIP_TRAF2
k34, IKK --> IKKa
(k35,l35), 0 <--> IKK
(k36,l36), 0 <--> NFkB
(k37,l37), 0 <--> FLIP
(k38,l38), 0 <--> XIAP
(k39,l39), 0 <--> A20
k40, IKKa --> 0
k41, IkBa_NFkB --> 0
k42, NFkB_N --> 0
k43, IkBa_mRNA --> 0
k44, IkBa --> 0
k45, IkBa_N --> 0
k46, IkBa_NFkB_N --> 0
k47, PIkBa --> 0
k48, A20_mRNA --> 0
k49, XIAP_mRNA --> 0
k50, FLIP_mRNA --> 0
k51, IKK --> IKKa
k52, IKKa --> IKK
k53, TNFRC1 --> TNF_TNFR_TRADD + TRAF2
k54, IkBa + NFkB --> IkBa_NFkB
k55, IkBa_NFkB --> NFkB + PIkBa
k56, NFkB --> NFkB_N
k57, 0 --> IkBa_mRNA
k58, 0 --> IkBa
(k59,l59), IkBa <--> IkBa_N
k60, IkBa_N + NFkB_N --> IkBa_NFkB_N
k61, IkBa_NFkB_N --> IkBa_NFkB
k62, 0 --> A20_mRNA
k63, 0 --> A20
k64, 0 --> XIAP_mRNA
k65, 0 --> XIAP
k66, 0 --> FLIP_mRNA
k67, 0 --> FLIP
(k68,l68), 0 <--> pCasp8
(k69,l69), 0 <--> pCasp3
(k70,l70), 0 <--> pCasp6
k71, Casp8 --> 0
k72, Casp3 --> 0
k73, Casp6 --> 0
k74, XIAP_Casp3 --> 0
(k75,l75), 0 <--> BAR
k76, BAR_Casp8 --> 0
(k77,l77), PARP <--> 0
k78, cPARP --> 0
k79, pCasp3 --> Casp3
k80, pCasp6 --> Casp6
k81, pCasp8 --> Casp8
(k82,l82), Casp3 + XIAP <--> XIAP_Casp3
k83, XIAP --> 0
k84, XIAP_Casp3 --> XIAP
k85, RIP --> 0
k86, FLIP --> 0
k87, PARP --> cPARP
(k88,l88), BAR + Casp8 <--> BAR_Casp8
end


#Set up the data
S = matrix(ZZ, netstoichmat(rn)) #stoichiometric matrix
Sminus = matrix(ZZ, substoichmat(rn)) #reactant matrix
T, _, R = get_symbolic_matrix(transpose(Sminus)) #T is the polynomial ring where the symbols in Rt live in
C = nullspace(S)[2]
@time ccf_res = ccf(R,C,columnLM = true)

@time for i in eachindex(ccf_res.Cblocks)
    diagonal_block = vcat(matrix(T,ccf_res.P[ccf_res.RblocksQ[i],ccf_res.Cblocks[i]]),
                          transpose(R)[ccf_res.RblocksT[i],ccf_res.Cblocks[i]])
    D = det(diagonal_block)
    println(length(ccf_res.Cblocks[i]),", ",length(vars(D)),", ", length(monomials(D)))
end

#did not finish in 30 minutes
#A = hcat(matrix(T,C),R)
#det(A)

#Column CCF
column_block_order(ccf_res)
column_block_order_pairs(ccf_res)
block_poset(ccf_res)
