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
