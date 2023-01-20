#we also add some function for 2d, so the name is old
using LinearAlgebra
export convert_from_spacial_and_spin_idx, get_periodic_idx,gene_H_k, gene_H_local, gene_ring_hopping,gene_2d_hopping, proj_K, proj_local, gene_ni, gene_nk, gene_n0σ, map_to_2d_idx

# length(basis)                   => 36
# size(ops)  => 8x8
"""
for the one band case, 
spin=0,1, (i.e up and down)
we store the spin up 1..L first then L+1, ..., 2L for spin down
convert_from_spacial_and_spin_idx(1,1,4)
convert_from_spacial_and_spin_idx(1,0,4)
"""
function convert_from_spacial_and_spin_idx(i,spin,L)
    i+spin*L
end

"""
get_periodic_idx(1,4)
get_periodic_idx(4,4)
get_periodic_idx(5,4)
get_periodic_idx(0,4)
get_periodic_idx(-1,4)
"""
function get_periodic_idx(idx,L)
    idx=idx%L
    if(idx<=0)
        idx+=L
    end
    idx
end

"""
this create the kinetic Hamiltonian with nearest hopping terms.
Hk=gene_H_k(L,ops)
"""
function gene_H_k(L,ops)
    H_k=zero(ops[1,1])
    for i in 1:L
        j=get_periodic_idx(i+1,L)
        i_up=convert_from_spacial_and_spin_idx(i,0,L)
        j_up=convert_from_spacial_and_spin_idx(j,0,L)
        i_dn=convert_from_spacial_and_spin_idx(i,1,L)
        j_dn=convert_from_spacial_and_spin_idx(j,1,L)
        H_k+=ops[i_up,j_up]+ops[i_dn,j_dn]+ops[j_up,i_up]+ops[j_dn,i_dn] 
   end
    H_k
end

"""
the double occupancy terms
"""
function gene_H_local(L,ops)
    H_local=zero(ops[1,1])
    for i in 1:L
        i_up=convert_from_spacial_and_spin_idx(i,0,L)
        i_dn=convert_from_spacial_and_spin_idx(i,1,L)
        H_local+=ops[i_up,i_up]*ops[i_dn,i_dn]
    end
    H_local
end

"""
create the non-interacting Hamiltonian, used to create the real k-point basis
"""
function gene_ring_hopping(L)
    result=zeros(L,L)
    for i in 1:L
        j=get_periodic_idx(i+1,L)
        result[i,j]=result[j,i]=1
    end
    result
end
"""
generate a hopping for 2d lattice
Lx=2
Ly=2
H_non_int=gene_2d_hopping(Lx,Ly)
ϵk,vk_trans=eigen(H_non_int)
H_non_int_1=gene_ring_hopping(Lx*Ly)
ϵk_1,vk_trans_1=eigen(H_non_int_1)

"""
function gene_2d_hopping(Lx,Ly)
    L=Lx*Ly
    result=zeros(L,L)
    for x_idx in 1:Lx
        for y_idx in 1:Ly
            idx=map_to_2d_idx(x_idx,y_idx,Lx,Ly)
            idx_x_plus_1=map_to_2d_idx(get_periodic_idx(x_idx+1,Lx),y_idx,Lx,Ly)
            idx_y_plus_1=map_to_2d_idx(x_idx,get_periodic_idx(y_idx+1,Ly),Lx,Ly)
            result[idx,idx_x_plus_1]=result[idx_x_plus_1,idx]=1
            result[idx,idx_y_plus_1]=result[idx_y_plus_1,idx]=1
        end
    end
    result
end

"""
conviert (x_idx,y_idx) -> idx
map_to_2d_idx(2,2,2,2)
"""
function map_to_2d_idx(x_idx,y_idx,Lx,Ly)
    (x_idx-1)*Lx+y_idx
end


# first, we build the basis and ops
# test1=Hk-(-2.0*nk_up[1]+2.0*nk_up[4]-2.0*nk_dn[1]+2.0*nk_dn[4])
# sum(abs.(test1))
# creat projector using λ_k
"""
for individual term
λ=0.5
n=nk[1]
"""
function proj_K(λ::Number,n)
    id=Matrix{Float64}(I,size(n))
    (1-λ)*id+(2*λ-1)*n
end  
  
"""
λk_σ=[1.0,0.5,0.5,0.0]
nk_σ=nk_up
"""
function proj_K(λk_σ::Array,nk_σ)
    N_orb=length(λk_σ)
    proj=Matrix{Float64}(I,size(nk_σ[1]))
    for k in 1:N_orb
        proj=proj*proj_K(λk_σ[k],nk_σ[k])
    end
    proj
end

"""
local projector for a given site, 
n_up=ni_up[1]
n_dn=ni_dn[1]
"""
function proj_local(u::Number,n_up,n_dn)
    id=Matrix{Float64}(I,size(n_up))
    p_up=id-n_up
    p_dn=id-n_dn
    X0=p_up*p_dn
    X2=n_up*n_dn
    Xup=n_up*p_dn
    Xdn=n_dn*p_up
    # sum(abs.(X0+Xup+Xdn+X2-id)) =>0  checked
    (1+u/4)*(X0+X2)+(1-u/4)*(Xup+Xdn)
end

"""
we pass an array of u 
u=[0.0 for _ in 1:4]
"""
function proj_local(u::Array,ni_up,ni_dn)
    N_orb=length(u)
    proj=Matrix{Float64}(I,size(ni_up[1]))
    for i in 1:N_orb
        proj=proj*proj_local(u[i],ni_up[i],ni_dn[i])
    end
    proj
end


"""
generate real space density, i.e, just the diagonal part of blocks
"""
function gene_ni(ops,spin,N_orb)
    ni=Array{Any,1}(undef,N_orb)
    for i in 1:N_orb
        i_spin=convert_from_spacial_and_spin_idx(i,spin,N_orb)
        ni[i]=ops[i_spin,i_spin]
    end
    ni
end



"""
compute the nk
k=1
# eigenvalues and eigenvector for non-interacting Hamiltonian
# vk[i,:] , vector for ith diagonal vector in real space basis, consistent with the formulas
"""
function gene_nk(vk,ops,spin)
    N_orb=size(vk)[1]
    nk=Array{Any,1}(undef,N_orb)
    for k in 1:N_orb
        v=vk[k,:]
        nk[k]=zero(ops[1,1])
        for i in 1:N_orb
            for j in 1:N_orb
                i_spin=convert_from_spacial_and_spin_idx(i,spin,N_orb)
                j_spin=convert_from_spacial_and_spin_idx(j,spin,N_orb)
                nk[k]+=v[i]*v[j]*ops[i_spin,j_spin]
            end
        end
    end
    nk
end

"""
λk_σ=λk_up
vi=vk_trans
"""
function gene_n0σ(λk_σ,vi)
    N_orb=size(vi)[1]
    n0σ=zeros(N_orb,N_orb)
    for i in 1:N_orb
        for j in 1:N_orb
            for k in 1:N_orb
                n0σ[i,j]+=vi[i,k]*vi[j,k]*λk_σ[k]
            end
        end
    end
    n0σ
end
