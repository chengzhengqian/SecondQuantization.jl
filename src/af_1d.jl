# auxilliary field
using Statistics
L=4
# seed1=1
# seed2=2
# g1=digits(seed1,base=2,pad=L)
# g2=digits(seed2,base=2,pad=L)

"""
Γ=g1
i=1
spin=0
u=-1
expv_from_config(Γ,spin,u)
"""
function expv_from_config(Γ,spin,u)
    expv_=(-u+4*sqrt(-u)+4)/(u+4)
    expvs_=[expv_,1/expv_]
    N_orb=length(Γ)
    result=zeros(N_orb,N_orb)
    for i in 1:N_orb
        idx=(Γ[i]+spin)%2
        result[i,i]=expvs_[idx+1]
    end
    result
end

# we first compute for a given spin
u=-0.0
spin=0
# for a given spin
n0σ=n0_up
g1=digits(seed1,base=2,pad=L)
g2=digits(seed2,base=2,pad=L)

function compute_n_from_expv(expv)
    inv(id+transpose(inv(expv)))
end

"""
to check the formulas for general case
expv1=rand(3,3)
expv2=rand(3,3)
expv0=rand(3,3)
z_test1,n_test1=cal_z_and_n_from_expv(expv1,expv2,expv0)
n0=compute_n_from_expv(expv0)
z_test2,n_test2=cal_z_and_n_from_n0(expv1,expv2,n0)
z_test1-z_test2
sum(abs.(n_test1-n_test2))
"""
function cal_z_and_n_from_expv(expv1,expv2,expv0)
    id=Matrix{Float64}(I,size(expv1))
    expv=expv1*expv0*expv2
    z1=det(id+expv)
    z0=det(id+expv0)
    z=z1/z0
    n=compute_n_from_expv(expv)
    return z,n
end

"""
check for the general formulas
n0=compute_n_from_expv(expv0)
"""
function cal_z_and_n_from_n0(expv1,expv2,n0)
    id=Matrix{Float64}(I,size(expv1))
    expv21=expv2*expv1
    z=det(id+(expv21-id)*transpose(n0))
    expv1t=transpose(expv1)
    expv2t=transpose(expv2)
    n=expv2t*inv(n0*expv1t*expv2t+(id-n0))*n0*expv1t
    return z,n
end


"""
z_up,n_up=cal_z_and_n(g1,g2,0,u,n0_up)
z_dn,n_dn=cal_z_and_n(g1,g2,1,u,n0_dn)
"""
function cal_z_and_n(g1,g2,spin,u,n0σ)
    expv1=expv_from_config(g1,spin,u)
    expv2=expv_from_config(g2,spin,u)
    # expv21=expv2*expv1
    # id=Matrix{Float64}(I,size(expv1))
    # zΓ=det(id+(expv21-id)*transpose(n0σ))
    # expv1T=transpose(expv1)
    # expv2T=transpose(expv2)
    # nΓ=expv2T*inv(n0σ*expv1T*expv2T+(id-n0σ))*n0σ*expv1T
    zΓ,nΓ=cal_z_and_n_from_n0(expv1,expv2,n0σ)
end



function sum_af_exact(u,L,n0_up,n0_dn)
    Z_total=0
    n_up_bare=zeros(L,L)
    n_dn_bare=zeros(L,L)
    dval_bare=zeros(L)
    for seed1 in 0:(2^L-1)
        for seed2 in 0:(2^L-1)
            g1=digits(seed1,base=2,pad=L)
            g2=digits(seed2,base=2,pad=L)
            # 0, spin up , 1, spin dn
            z_up,n_up=cal_z_and_n(g1,g2,0,u,n0_up)
            z_dn,n_dn=cal_z_and_n(g1,g2,1,u,n0_dn)
            dval=[n_up[i,i]*n_dn[i,i] for i in 1:L]
            z=z_up*z_dn
            Z_total+=z
            n_up_bare+=n_up*z
            n_dn_bare+=n_dn*z
            dval_bare+=dval*z
        end
    end
    n_up_bare/Z_total, n_dn_bare/Z_total,dval_bare/Z_total
end

function gene_data_af_exact(u)
    n_up,n_dn,dval=sum_af_exact(u,L,n0_up,n0_dn)
    u,mean(dval),tr(H_non_int*(n_up+n_dn))/L
end

L=4
N_up=2
N_dn=2
filename_tag="L_$(L)_N_up_$(N_up)_N_dn_$(N_dn)"
H_non_int=gene_ring_hopping(L)
ϵk,vk_trans=eigen(H_non_int)
λk_up=[1.0,0.5,0.5,0.0]
λk_dn=[1.0,0.5,0.5,0.0]
n0_up=gene_n0σ(λk_up,vk_trans)
n0_dn=gene_n0σ(λk_dn,vk_trans)

us=-3.5:0.1:0.0
data=[gene_data_af_exact(u_) for u_ in us]
saveData(data,"$(data_dir)/$(filename_tag)_u_d_k_af_exact.dat")


"""
we first verify the two formulas to compute z and n
"""

function cal_λ_from_ϵ(ϵ,μ)
    if (ϵ>μ+1E-5)
        0.0
    elseif(ϵ<μ-1E-5)
        1.0
    else
        0.5
    end
end


L=100
N_up=50
N_dn=50
filename_tag="L_$(L)_N_up_$(N_up)_N_dn_$(N_dn)"
H_non_int=gene_ring_hopping(L)
ϵk,vk_trans=eigen(H_non_int)
λk_up=(ϵ->cal_λ_from_ϵ(ϵ,0.0)).(ϵk)
λk_dn=(ϵ->cal_λ_from_ϵ(ϵ,0.0)).(ϵk)
n0_up=gene_n0σ(λk_up,vk_trans)
n0_dn=gene_n0σ(λk_dn,vk_trans)
# n_up_test,n_dn_test,dval_test=sum_af_monte_carlo(-0.5,L,n0_up,n0_dn,1000000)
# n_up_test1,n_dn_test1,dval_test1=sum_af_monte_carlo(-0.5,L,n0_up,n0_dn,1000)
# mean(dval_test)
# mean(dval_test1)
us=-3.0:0.1:0.0
N_iter=100000
tag="diagonal"
N_iter=10000
N_iter=100000
tag="diagonal_sub"
data=[gene_data_af_monte_carlo(u_,tag,N_iter) for u_ in us]
# saveData(data,"$(data_dir)/$(filename_tag)_u_d_k_af_monte_carlo.dat")
# saveData(data,"$(data_dir)/$(filename_tag)_u_d_k_af_monte_carlo_N_iter_$(N_iter).dat")
saveData(data,"$(data_dir)/$(filename_tag)_u_d_k_af_monte_carlo_N_iter_$(N_iter)_tag_$(tag).dat")

g_anti=Array{Int64,1}(undef,L)
g_anti[1:2:99].=0
g_anti[2:2:100].=1
u=-2.0
g_random=rand([0,1],L)
z_up_new,n_up_new=cal_z_and_n(g_anti,g_anti,0,u,n0_up)
z_dn_new,n_dn_new=cal_z_and_n(g_anti,g_anti,1,u,n0_dn)
dval=[n_up_new[i,i]*n_dn_new[i,i] for i in 1:L]
mean(dval)
z_up_random,n_up_random=cal_z_and_n(g_random,g_random,0,u,n0_up)
z_dn_random,n_dn_random=cal_z_and_n(g_random,g_random,1,u,n0_dn)
dval_random=[n_up_random[i,i]*n_dn_random[i,i] for i in 1:L]
mean(dval_random)
den_up=[n_up_random[i,i] for i in 1:L]
dval_random[end-10]
g_random[(end-20):(end-5)]
den_up[(end-20):(end-5)]



function gene_data_af_monte_carlo(u,tag,N_iter=10000)
    n_up,n_dn,dval=sum_af_monte_carlo(u,L,n0_up,n0_dn,N_iter,tag)
    u,mean(dval),tr(H_non_int*(n_up+n_dn))/L
end


# now, we consider the Monte Carlo sampling
"""
N_iter=10000
seed=1
n_up,n_dn,dval=sum_af_exact(-3.0,L,n0_up,n0_dn)
n_up_test,n_dn_test,dval_test=sum_af_monte_carlo(-2.0,L,n0_up,n0_dn,10000)
# try some large system
mean(dval_test)
tr(H_non_int*(n_up+n_dn))/L
tr(H_non_int*(n_up_test+n_dn_test))/L
u=-1.0
tag="random"
tag="diag"
see generate_new_config
"""
function sum_af_monte_carlo(u,L,n0_up,n0_dn,N_iter,tag)
    N_total=0
    Z_total=0
    N_new=0
    n_up_bare=zeros(L,L)
    n_dn_bare=zeros(L,L)
    dval_bare=zeros(L)
    g_anti=Array{Int64,1}(undef,L)
    g_anti[1:2:99].=0
    g_anti[2:2:100].=1
    g1,g2=generate_new_config(g_anti,g_anti,tag)
    z_up,n_up=cal_z_and_n(g1,g2,0,u,n0_up)
    z_dn,n_dn=cal_z_and_n(g1,g2,1,u,n0_dn)
    z=z_up*z_dn
    p=abs(z)
    # rand([1,2],10,10)
    for i in 1:N_iter
        if(i%10==0)
            print("$(i) iter\n")
        end
        g1_new,g2_new=generate_new_config(g1,g2,tag)
        z_up_new,n_up_new=cal_z_and_n(g1_new,g2_new,0,u,n0_up)
        z_dn_new,n_dn_new=cal_z_and_n(g1_new,g2_new,1,u,n0_dn)
        z_new=z_up_new*z_dn_new
        p_new=abs(z_new)
        α=rand(1)[1]
        if(p_new/p>=α)
            # take the new sample
            z=z_new
            p=p_new
            n_up=n_up_new
            n_dn=n_dn_new
            N_new+=1
        end
        dval=[n_up[i,i]*n_dn[i,i] for i in 1:L]
        # as we sample according to p
        N_total+=1
        Z_total+=sign(z)
        n_up_bare+=n_up*sign(z)
        n_dn_bare+=n_dn*sign(z)
        dval_bare+=dval*sign(z)
    end
    print("N_total $(N_total), Z_total: $(Z_total), N_new/N_total:$(N_new/N_total)\n")
    n_up_bare/Z_total, n_dn_bare/Z_total,dval_bare/Z_total
end

"""
we randomly choose g1 or g2, and randomly flip one of the bit
"""
function generate_new_config(g1,g2,tag)
    if(tag=="random")
        g1=rand([0,1],L)
        g2=rand([0,1],L)
        return (g1,g2)
    elseif(tag=="diagonal")
        g1=rand([0,1],L)
        return (g1,g1)
    elseif(tag=="diagonal_sub")
        idx=rand(1:L,1)[1]
        g1[idx]=rand([0,1],1)[1]
        return (g1,g1)
    else
        error("tag $(tag) is not supported!\n")
    end
end
