]activate "/home/chengzhengqian/share_workspace/czq_julia_package/SecondQuantization"

using Revise
using SecondQuantization

# create state
s1=QuantumState("0001")
s2=QuantumState("1010")

op1=op_Num(1)
op1(s2)

getNum(s2,5)
# 
# SecondQuantization.gene_basis_and_op
# gene_basis_and_op
# test the repr
L=4
N_up=N_dn=2
basis,ops=gene_basis_and_op(L,N_up,N_dn)
# size(ops[1])
# length(basis.IdxToState)
N_orbitals=2*L
# SecondQuantization.repr(op_Create(1),basis)
# check the op_Create, we need a full fock space basis
basis,ops=gene_basis_and_op(L)
size(ops)
size(ops[1])
cr_ops=[SecondQuantization.repr(op_Create(i),basis) for i in 1:L]
an_ops=[SecondQuantization.repr(op_Annihilate(i),basis) for i in 1:L]
ops_alt=Array{Any,2}(undef,L,L)
errors=zeros(L,L)
for i in 1:L
    for j in 1:L
        ops_alt[i,j]=cr_ops[i]*an_ops[j]
        errors[i,j]=sum(abs.(ops_alt[i,j]-ops[i,j]))
    end
end

i=j=1
i=1;j=2
cr_ops[i]*an_ops[j]+an_ops[j]*cr_ops[i]




