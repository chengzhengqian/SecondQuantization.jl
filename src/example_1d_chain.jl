using LinearAlgebra
# example of using 1d chain
# benchmark for auxillary field

# we first test gutzwiler, N=2

# we first generate the N=2 exact solution
using SecondQuantization
function gene_data_exact(u_)
    u=[u_ for _ in 1:L]
    proj1=proj_local(u,ni_up,ni_dn)
    ρ=proj1*ρ0*proj1
    dval=expt(Hloc,ρ)/L
    kval=expt(Hk,ρ)/L
    dval,kval
end

L=4
N_up=2
N_dn=2
filename_tag="L_$(L)_N_up_$(N_up)_N_dn_$(N_dn)"
basis,ops=gene_basis_and_op(L,N_up,N_dn)
Hk=gene_H_k(L,ops)
Hloc=gene_H_local(L,ops)
# we first study the N=2 G-type
H_non_int=gene_ring_hopping(L)
ϵk,vk_trans=eigen(H_non_int)
vk=transpose(vk_trans)
nk_up=gene_nk(vk,ops,0)
nk_dn=gene_nk(vk,ops,1)
ni_up=gene_ni(ops,0,L)
ni_dn=gene_ni(ops,1,L)
# this part are fixed
λk_up=[1.0,0.5,0.5,0.0]
λk_dn=[1.0,0.5,0.5,0.0]
# single particle density matrix from λk_σ
n0_up=gene_n0σ(λk_up,vk_trans)
n0_dn=gene_n0σ(λk_dn,vk_trans)
projK_up=proj_K(λk_up,nk_up)
projK_dn=proj_K(λk_dn,nk_dn)
ρ0=projK_up*projK_dn



data_dir="./1d_chain_data"
mkdir(data_dir)
us=-4.0:0.01:0.0
data=[(u_,gene_data_exact(u_)...) for u_ in us]
using CZQUtils
saveData(data,"$(data_dir)/$(filename_tag)_u_d_k.dat")


L=4
N_up=1
N_dn=1
filename_tag="L_$(L)_N_up_$(N_up)_N_dn_$(N_dn)"
basis,ops=gene_basis_and_op(L,N_up,N_dn)
Hk=gene_H_k(L,ops)
Hloc=gene_H_local(L,ops)
# we first study the N=2 G-type
H_non_int=gene_ring_hopping(L)
ϵk,vk_trans=eigen(H_non_int)
vk=transpose(vk_trans)
nk_up=gene_nk(vk,ops,0)
nk_dn=gene_nk(vk,ops,1)
ni_up=gene_ni(ops,0,L)
ni_dn=gene_ni(ops,1,L)
# this part are fixed
λk_up=[1.0,0.0,0.0,0.0]
λk_dn=[1.0,0.0,0.0,0.0]
# single particle density matrix from λk_σ
n0_up=gene_n0σ(λk_up,vk_trans)
n0_dn=gene_n0σ(λk_dn,vk_trans)
projK_up=proj_K(λk_up,nk_up)
projK_dn=proj_K(λk_dn,nk_dn)
ρ0=projK_up*projK_dn



data_dir="./1d_chain_data"
mkdir(data_dir)
us=-4.0:0.01:0.0
data=[(u_,gene_data_exact(u_)...) for u_ in us]
using CZQUtils
saveData(data,"$(data_dir)/$(filename_tag)_u_d_k.dat")



