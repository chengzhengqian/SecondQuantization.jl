export gene_basis_and_op, exp_c, expt, compute_density_matrix,dot_c
# move these useful function into a utility file
# these create the operator represetation
"""
N is the number of spin orbitals,
ops is matrix the single particle operator
"""
function gene_basis_and_op(N)
    basis=gene_basis(N)
    ops=gene_ops(basis,N)
    basis, ops
end

function gene_basis_and_op(N,N_number)
    basis=gene_basis(N,N_number)
    ops=gene_ops(basis,N)
    basis, ops
end
"""
N is the number of orbital, and 2*N is the total number of spin-orbitals
"""
function gene_basis_and_op(N,N_up,N_dn)
    basis=gene_basis(N,N_up,N_dn)
    ops=gene_ops(basis,2*N)
    basis, ops
end

function gene_ops(basis,N)
    ops=Array{Any,2}(undef,N,N)
    for i in 1:N
        for j in 1:N
            if(i!=j)
                ops[i,j]=repr(op_Hopping(j,i),basis) #  a_d_i*aj
            else
                ops[i,j]=repr(op_Num(i),basis)
            end
        end
    end
    ops
end
# compute the exp(cij*nij)
function exp_c(c,ops)
    N=size(c)[1]
    result=zeros(2^N,2^N)
    for i in 1:N
        for j in 1:N
            result+=c[i,j]*ops[i,j]
        end
    end
    exp(result)
end

function dot_c(c,ops)
    N=size(c)[1]
    result=zeros(2^N,2^N)
    for i in 1:N
        for j in 1:N
            result+=c[i,j]*ops[i,j]
        end
    end
    result
    # exp(result)
end

"""
compute expectation values
"""
function expt(A,rho)
    tr(A*rho)/tr(rho)
end

function compute_density_matrix(rho,ops)
    N_=size(ops)[1]
    C=tr(rho)
    result=zeros(N_,N_)
    for i in 1:N_
        for j in 1:N_
            result[i,j]=tr(rho*ops[i,j])
        end
    end
    result/C
end

function compute_density_matrix_half(rho,ops)
    N_=toInt(size(ops)[1]/2)
    C=tr(rho)
    result=zeros(N_,N_)
    for i in 1:N_
        for j in 1:N_
            result[i,j]=tr(rho*ops[i,j])
        end
    end
    result/C
end

function compute_density_density_correlation(rho,ops)
    N_=size(ops)[1]
    C=tr(rho)
    result=zeros(N_,N_)
    for i in 1:N_
        for j in 1:N_
            result[i,j]=tr(rho*ops[i,i]*ops[j,j])
        end
    end
    result/C
end
