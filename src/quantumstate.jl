# provide a simple interface of the package
using Base
using CZQUtils

export op, gene_basis, QuantumState, op_Num, op_Hopping, op_NN, op_p, gene_basis, getNum, op_Create,op_Annihilate, repr
    


mutable struct  QuantumState
    state::UInt64
    size::Unsigned
    coef
    QuantumState(state,size)=new(state,size,1.0)
    QuantumState(state,size,coef)=new(state,size,coef)
end


function QuantumState(str::String)
    size=length(str)
    new_state=QuantumState(0x0,size,1.0)
    for i in 1:size
        if(str[i]=='0')
            setNum!(new_state,i,0)
        else
            setNum!(new_state,i,1)
        end
    end
    new_state
end
using Base

Base.hash(x::QuantumState, h::UInt)=hash(x.state,hash(x.size,h))

Base.isequal(x::QuantumState, y::QuantumState)=Base.isequal(x.state,y.state)&&Base.isequal(x.size,y.size)

function Base.show(io::IO, x::QuantumState)
    print("( ");print(x.coef);print(" |")
    for i in 1:x.size
        num=getNum(x,i)
        if(num==0)
            print("0")
        else
            print("1")
        end
    end
    print(">)")
end

# x=QuantumState(8,5,1.0)

# using Base
# using JuliaCommon
# 1&3
"""
we should raise exception, in case we may wrongly think the idx is varid
"""
function getSafeIdx(x::QuantumState, idx)
    if(idx<1)
        error("index $idx <1 \n")
        idx=1
    end
    if(idx>x.size)
        error("index $idx > $(x.size) \n")
        idx=x.size
    end
    idx
end

function getNum(x::QuantumState, idx)
    idx=getSafeIdx(x,idx)
    (x.state&(UInt64(0x1)<<(idx-1)))>>(idx-1)
end
# total number
function getNum(x::QuantumState)
   sum([getNum(x,idx) for idx in 1:x.size])
end

# getNum(QuantumState("110111"))
# get the particle number in between idx1 (exclude) idx2 (exclue)

# min(1,2)
"""
compute the number between idx1 and idx2, (idx1<=idx2) after sort, not including idx2 and idx1
s=QuantumState("1010")
getNum(s,0,1) ->0
getNum(s,0,2) ->1
"""
function getNum(x::QuantumState, idx1,idx2)
    # idx1=getSafeIdx(x,idx1);idx2=getSafeIdx(x,idx2)
    min_idx=min(idx1,idx2);max_idx=max(idx1,idx2)
    nums=[getNum(x,idx) for idx in(min_idx+1):(max_idx-1)]
    if(length(nums)==0)
        return 0
    end
    sum(nums)
end


function setNum(state::Number,idx,num)
    if(num==0)
        new_state=state&(~(UInt64(0x1)<<(idx-1)))
    else
        new_state=state|(UInt64(0x1)<<(idx-1))
    end
    new_state
end

function setNum(x::QuantumState, idx,num)
    idx=getSafeIdx(x,idx)
    QuantumState(setNum(x.state,idx,num),x.size,x.coef)
end

function mulCoef(x::QuantumState,factor)
    QuantumState(x.state,x.size,x.coef*factor)
end

function mulCoef!(x::QuantumState,factor)
    x.coef=x.coef*factor
end

function setNum!(x::QuantumState, idx,num)
    idx=getSafeIdx(x,idx)
    x.state=setNum(x.state,idx,num)
    x
end

# x
# setNum(x,2,1)
struct Op
    op::Function
end

(f::Op)(x)=f.op(x)

"""
a_dag_idx*a_idx
"""
function op_Num(idx)
    op=(x::QuantumState)->mulCoef(x,getNum(x,idx))
    Op(op)
end

"""
we also add the op for creation and annihilation operator
x=QuantumState("1011")
idx=2
op1=op_Create(2)
op1(x)
x=QuantumState("0011")
op1(x)
x=QuantumState("0111")
op1(x)
"""
function op_Create(idx)
    function op(x::QuantumState)
        num_before_idx=Float64(getNum(x,0,idx))
        num_at_idx=Float64(getNum(x,idx))
        factor=(1-num_at_idx)*(-1)^num_before_idx
        new_state=mulCoef(x,factor)
        setNum!(new_state,idx,1)
        new_state
    end
    Op(op)
end

"""
Annihilate the particle at site idx
op2=op_Annihilate(2)
s
op2(op1(op2(op1(s))))
"""
function op_Annihilate(idx)
    function op(x::QuantumState)
        num_before_idx=Float64(getNum(x,0,idx))
        num_at_idx=Float64(getNum(x,idx))
        factor=num_at_idx*(-1)^num_before_idx
        new_state=mulCoef(x,factor)
        setNum!(new_state,idx,0)
        new_state
    end
    Op(op)
end


"""
a_dag_idx2*a_idx1
"""
function op_Hopping(idx1,idx2)
    function op(x::QuantumState)
        num1=Float64(getNum(x,idx1))
        num2=Float64(getNum(x,idx2))
        num_between=getNum(x,idx1,idx2)
        factor=num1*(1-num2)*(-1)^num_between
        new_state=mulCoef(x,factor)
        if(factor!=0.0)
            setNum!(new_state,idx1,0)
            setNum!(new_state,idx2,1)
        end
        new_state
    end
    Op(op)
end


# n1*n2*n3..
"""
n_idx1*n_idx2*..
where idx_i is in the array of idx
"""
function op_Num(idx::AbstractArray)
    function op(x::QuantumState)
        new_state=mulCoef(x,getNum(x,idx[1]))
        for idx_ in idx[2:end]
            mulCoef!(new_state,getNum(x,idx_))
        end
        new_state
    end
    Op(op)
end
# (hat{n}-nc)

# function op_NN(idx::AbstractArray,nc::Number)
#     function op(x::QuantumState)
#         new_state=mulCoef(x,getNum(x,idx[1])-nc)
#         for idx_ in idx[2:end]
#             mulCoef!(new_state,getNum(x,idx_)-nc)
#         end
#         new_state
#     end
#     Op(op)
# end

function op_NN(idx::AbstractArray,nc::Number)
    function op(x::QuantumState)
        new_state=mulCoef(x,getNum(x,idx[1])-nc)
        for idx_ in idx[2:end]
            mulCoef!(new_state,getNum(x,idx_)-nc)
        end
        new_state
    end
    Op(op)
end
# create density matrix from local reduced density matrix
# assuming e+up+dn+d=1

function op_ρ(e,up,dn,d)
    function op(x::QuantumState)
        N_lattice=Int(x.size/2)
        factor=1.0
        for i in 1:N_lattice
            n_up=getNum(x,i)
            n_dn=getNum(x,i+N_lattice)
            if(n_up==1&&n_dn==1)
                factor*=d
            elseif(n_up==1&&n_dn==0)
                factor*=up
            elseif(n_up==0&&n_dn==1)
                factor*=dn
            else
                factor*=e
            end
        end
        mulCoef(x,factor)
    end
    Op(op)
end

# b_test1=gene_basis(4)
# Int(get(b_test1,1).size/2)
# sum(repr(op_ρ(0.1,0.4,0.4,0.1),b_test1))

mutable struct QuantumBasis
    IdxToState::Dict{Integer,QuantumState}
    StateToIdx::Dict{QuantumState,Integer}
end

function addToBasis!(basis,state,index)
    basis.StateToIdx[state]=index
    basis.IdxToState[index]=state
end
function get(basis::QuantumBasis,idx::Number)
    basis.IdxToState[idx]
end
function get(basis::QuantumBasis,state::QuantumState)
    basis.StateToIdx[state]
end
function has(basis::QuantumBasis,idx::Number)
    haskey(basis.IdxToState,idx)
end
function has(basis::QuantumBasis,state::QuantumState)
    haskey(basis.StateToIdx,state)
end
Base.length(basis::QuantumBasis)=length(basis.IdxToState)
# general style
function gene_basis(N_lattice)
    basis=QuantumBasis(Dict(),Dict())
    index=1
    for i in 0:(2^N_lattice-1)
        state=QuantumState(i,N_lattice)
        addToBasis!(basis,state,index)
        index+=1
    end
    basis
end
    
# spinless first
function gene_basis(N_lattice,N_number)
    basis=QuantumBasis(Dict(),Dict())
    index=1
    for i in 0:(2^N_lattice-1)
        state=QuantumState(i,N_lattice)
        if(getNum(state)==N_number)
            addToBasis!(basis,state,index)
            index+=1
        end
    end
    basis
end
# spin one band
function gene_basis(N_lattice,N_up,N_dn)
    basis=QuantumBasis(Dict(),Dict())
    index=1
    for i in 0:(2^(N_lattice*2)-1)
        state=QuantumState(i,N_lattice*2)
        if(getNum(state,0,N_lattice+1)==N_up
           && getNum(state,N_lattice,N_lattice*2+1)==N_dn
           )
            addToBasis!(basis,state,index)
            index+=1
        end
    end
    basis
end

# gene_basis(2,1,1)


# s=QuantumState(10,5)
# getNum(s,0,3)
# gene_basis(3,0,0)
# getNum(QuantumState("100001"),1,2)
# state=QuantumState(4,4)
# getNum(state,0,5)
# basis=gene_basis(4,2)
# n1=op_Num(1)
# get(basis,get(basis,2))
# length(basis)
# repr(n1,basis)
# repr(op_Hopping(1,4),basis)
# op(get(basis,1))
# h=repr(op_Hopping(1,4),basis)
# h=repr(op_Hopping(2,4),basis,h)
# s=get(basis,1)
# op=op_Hopping(1,4)
# op(s)
# s=QuantumState("1010")
# op(s)
"""
returns a zero matrix with size of basis
"""
function repr(basis::QuantumBasis)
    size=length(basis)
    zeros(size,size)
end

function repr(op::Op,basis::QuantumBasis;factor=1.0)
    result=repr(basis)
    for i in 1:length(basis)
        new_state=op(get(basis,i))
        if(has(basis,new_state))
            j=get(basis,new_state)
            result[j,i]=new_state.coef*factor
        else
            print("$(new_state) is no in basis")
        end
    end
    result
end

function repr(op::Op,basis::QuantumBasis,result;factor=1.0)
    for i in 1:length(basis)
        new_state=op(get(basis,i))
        if(has(basis,new_state))
            j=get(basis,new_state)
            result[j,i]+=new_state.coef*factor
        end
    end
    result
end



# basis1=gene_basis(10,5)
# length(basis1)
# basis
    
# s1=QuantumState("0001")
# s2=QuantumState("0010")
# a=Dict(s1=>1,s2=>2)
# h=op_Hopping(1,3)
# n=op_Num(6)
# n=op_Num([1,2])
# a[n(s1)]
# s=QuantumState("100110")
# s=QuantumState("111110")
# # n(s)
# n1=op_Num(3)
# n1(x)



