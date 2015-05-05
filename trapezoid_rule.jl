function trapezoid_rule(a,b,N)
    w=Array(Float64,(N+1,));
    x=linspace(a,b,N+1);
    w[1]=1.0;
    w[2:N]=2.0;
    w[N+1]=1.0;
    w[:]=w*((b-a)/(2*N));
    return (x,w);
end

function trapezoid_tensor_product(a,b,N)
    (x,w)=trapezoid_rule(a,b,N);
    y=copy(x);
    z=copy(x);
    W=Array(Float64,(N+1,N+1,N+1));
    for i = 1 : length(w)
        for j = 1 : length(w)
            for k = 1 : length(w)
                W[i,j,k]=w[i]*w[j]*w[k];
            end
        end
    end
    return (x,y,z,W);
end
