function sparse_reorder(A)
    (m,n)=size(A);
    U=Array(Array{Int64},(m,));
    for i=1:m
        U[i]=A.rowval[A.colptr[i]:A.colptr[i+1]-1];
    end
    p=sortperm(U,lt=lexless);
    return p;
end
