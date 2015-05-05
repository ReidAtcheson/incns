include("trapezoid_rule.jl");

 function ordertolen(order::Int64)
    #Three dimensional tensor product basis.
    expsz=3^order;

    #Using vector of three components.
    expsz=3*expsz;

    return expsz;
end

 function evaluate_tensor_basis(i::Int64,j::Int64,k::Int64,x::Float64,y::Float64,z::Float64)
    return exp(im*(i*x+j*y+k*z));
end

 function evaluate_tensor_basis_derivative(i::Int64,j::Int64,k::Int64,d::Int64,x::Float64,y::Float64,z::Float64)
    out=exp(im*(i*x+j*y+k*z));
    is=[i,j,k];
    out=im*is[d]*out;
    return out;
end
 function evaluate_tensor_basis_derivative(i::Int64,j::Int64,k::Int64,d1::Int64,d2::Int64,x::Float64,y::Float64,z::Float64)
    out=exp(im*(i*x+j*y+k*z));
    is=[i,j,k];
    out=(im*is[d2])*(im*is[d1])*out;
    return out;
end





 function evaluate_curl_basis_gradient(order::Int64,n::Int64,x::Float64,y::Float64,z::Float64)
    ncomps=3;
    m=div(n,ncomps);
    r=n%ncomps;
    (id1,id2,id3)=ind2sub((order,order,order),m);
    out=Array(Complex128,(ncomps,ncomps));
    if r==0
        dyx=evaluate_tensor_basis_derivative(id1,id2,id3,2,1,x,y,z);
        dzx=evaluate_tensor_basis_derivative(id1,id2,id3,3,1,x,y,z);
        dyy=evaluate_tensor_basis_derivative(id1,id2,id3,2,2,x,y,z);
        dzy=evaluate_tensor_basis_derivative(id1,id2,id3,3,2,x,y,z);
        dyz=evaluate_tensor_basis_derivative(id1,id2,id3,2,3,x,y,z);
        dzz=evaluate_tensor_basis_derivative(id1,id2,id3,3,3,x,y,z);
        out[:,1]=[0.0*im,dzx,-dyx];
        out[:,2]=[0.0*im,dzy,-dyy];
        out[:,3]=[0.0*im,dzy,-dyy];
    elseif r==1
        dxx=evaluate_tensor_basis_derivative(id1,id2,id3,1,1,x,y,z);
        dzx=evaluate_tensor_basis_derivative(id1,id2,id3,3,1,x,y,z);
        dxy=evaluate_tensor_basis_derivative(id1,id2,id3,1,2,x,y,z);
        dzy=evaluate_tensor_basis_derivative(id1,id2,id3,3,2,x,y,z);
        dxz=evaluate_tensor_basis_derivative(id1,id2,id3,1,3,x,y,z);
        dzz=evaluate_tensor_basis_derivative(id1,id2,id3,3,3,x,y,z);

        out[:,1]=[-dzx,0.0*im,dxx];
        out[:,2]=[-dzy,0.0*im,dxy];
        out[:,3]=[-dzz,0.0*im,dxz];
    else
        dxx=evaluate_tensor_basis_derivative(id1,id2,id3,1,1,x,y,z);
        dyx=evaluate_tensor_basis_derivative(id1,id2,id3,2,1,x,y,z);
        dxy=evaluate_tensor_basis_derivative(id1,id2,id3,1,2,x,y,z);
        dyy=evaluate_tensor_basis_derivative(id1,id2,id3,2,2,x,y,z);
        dxz=evaluate_tensor_basis_derivative(id1,id2,id3,1,3,x,y,z);
        dyz=evaluate_tensor_basis_derivative(id1,id2,id3,2,3,x,y,z);
        out[:,1]=[dyx,-dxx,0.0*im];
        out[:,2]=[dyy,-dxy,0.0*im];
        out[:,3]=[dyz,-dxz,0.0*im];
    end
    return out;
end



 function evaluate_curl_basis(order::Int64,n::Int64,x::Float64,y::Float64,z::Float64)
    ncomps=3;
    m=div(n,ncomps);
    r=n%ncomps;
    (id1,id2,id3)=ind2sub((order,order,order),m);
    out=(0.0*im,0.0*im,0.0*im);
    if r==0
        dy=evaluate_tensor_basis_derivative(id1,id2,id3,2,x,y,z);
        dz=evaluate_tensor_basis_derivative(id1,id2,id3,3,x,y,z);
        out=(0.0*im,dz,-dy);
    elseif r==1
        dx=evaluate_tensor_basis_derivative(id1,id2,id3,1,x,y,z);
        dz=evaluate_tensor_basis_derivative(id1,id2,id3,3,x,y,z);
        out=(-dz,0.0*im,dx);
    else
        dx=evaluate_tensor_basis_derivative(id1,id2,id3,1,x,y,z);
        dy=evaluate_tensor_basis_derivative(id1,id2,id3,2,x,y,z);
        out=(dy,-dx,0.0*im);
    end
    return out;
end

 function curl_basis_inner_product(order::Int64,m::Int64,n::Int64)
    ncomps=3;
    rm=m%ncomps;
    rn=n%ncomps;
    ordm=div(m,ncomps);
    ordn=div(n,ncomps);

    (id1,id2,id3)=ind2sub((order,order,order),ordm);
    idsm=[id1,id2,id3];
    (id1,id2,id3)=ind2sub((order,order,order),ordn);
    idsn=[id1,id2,id3];

    val=Float64((2*pi)*(2*pi)*(2*pi));

    rm=rm+1;
    rn=rn+1;
    if idsm==idsn
        if rm==rn
            if rm==1
                val=val*(idsn[3]*idsn[3]+idsn[2]*idsn[2])
            elseif rm==2
                val=val*(idsn[3]*idsn[3]+idsn[1]*idsn[1]);
            else 
                val=val*(idsn[2]*idsn[2]+idsn[1]*idsn[1]);
            end
        else
            val=-val*idsn[rn]*idsm[rm];
        end
    else
        val=0.0;
    end

    return val;
end

 function curl_basis_mass_matrix(order::Int64,quadlen::Int64)
    n=ordertolen(order);
    M=Array(Complex128,(n,n));
    (x,y,z,W)=trapezoid_tensor_product(-pi,pi,quadlen);
    V1=copy(W)*im;
    V2=copy(V1);
    V3=copy(V1);
    U1=copy(V1);
    U2=copy(V1);
    U3=copy(V1);
 
    for i = 1 : n
        for j = 1 : n
            #Evaluate basis at quadrature nodes.
            for u1 = 1 : length(x)
                for u2 = 1 : length(y)
                    for u3 = 1 : length(z)
                        (vx,vy,vz)=evaluate_curl_basis(order,i,x[u1],y[u2],z[u3]);
                        (ux,uy,uz)=evaluate_curl_basis(order,j,x[u1],y[u2],z[u3]);
                        V1[u1,u2,u3]=vx;
                        V2[u1,u2,u3]=vy;
                        V3[u1,u2,u3]=vz;

                        U1[u1,u2,u3]=ux;
                        U2[u1,u2,u3]=uy;
                        U3[u1,u2,u3]=uz;
     
                    end
                end
            end

            S=conj(V1).*U1+conj(V2).*U2+conj(V3).*U3;

            M[i,j]=sum(S[:].*W[:]);
        end
    end

    return M;
end

 function filter_dependent_modes(order::Int64,quadlen::Int64)
    M=curl_basis_mass_matrix(order,quadlen);
    r=rank(M);
    (Q,R,p)=qr(M,Val{true});
    return p[1:r];
end




 function filter_dependent_modes(order::Int64)
    M=full(curl_basis_mass_matrix(order));
    r=rank(M);
    (Q,R,p)=qr(M,Val{true});
    return p[1:r];
end

 function filter_dependent_modes_new(order::Int64)
    M=curl_basis_mass_matrix(order);
    (m,n)=size(M);
    S=IntSet([1:m;]);
    for s in S
      for t in S
	if t!=s
	  if M.rowval[M.colptr[t]:M.colptr[t+1]-1]==M.rowval[M.colptr[s]:M.colptr[s+1]-1]
	    v=(norm(M[:,s]'*M[:,t]))/(norm(M[:,s])*norm(M[:,t]));
	    if abs(v-1.0)<1e-3
	      print("eliminated pair=($s,$t)\n");
	      print("diff=($(s-t))\n");
	      print("rel diff=($((s-t)/n))\n");
	      delete!(S,t);
	    end
	  end
	end
      end
    end
    return [Int64(p) for p in S];
end


 function filter_dependent_modes_new(order::Int64,alpha)
    M=curl_basis_mass_matrix(order);
    (m,n)=size(M);
    S=IntSet([1:m;]);
    for s in S
      for t in S
	if t!=s
	  v=(norm(M[:,s]'*M[:,t]))/(norm(M[:,s])*norm(M[:,t]));
	  if abs(v-1.0)<alpha
	    print("eliminated pair=($s,$t)\n");
	    print("diff=($(s-t))\n");
	    print("rel diff=($((s-t)/n))\n");
	    delete!(S,t);
	  end
	end
      end
    end
    return [Int64(p) for p in S];
end






 function curl_basis_mass_matrix(p::Array{Int64},order::Int64,quadlen::Int64)
    n=length(p);
    M=Array(Complex128,(n,n));
    (x,y,z,W)=trapezoid_tensor_product(-pi,pi,quadlen);
    V1=copy(W)*im;
    V2=copy(V1);
    V3=copy(V1);
    U1=copy(V1);
    U2=copy(V1);
    U3=copy(V1);
 
    for i = 1 : n
        for j = 1 : n
            #Evaluate basis at quadrature nodes.
            for u1 = 1 : length(x)
                for u2 = 1 : length(y)
                    for u3 = 1 : length(z)
                        (vx,vy,vz)=evaluate_curl_basis(order,p[i],x[u1],y[u2],z[u3]);
                        (ux,uy,uz)=evaluate_curl_basis(order,p[j],x[u1],y[u2],z[u3]);
                        V1[u1,u2,u3]=vx;
                        V2[u1,u2,u3]=vy;
                        V3[u1,u2,u3]=vz;

                        U1[u1,u2,u3]=ux;
                        U2[u1,u2,u3]=uy;
                        U3[u1,u2,u3]=uz;
     
                    end
                end
            end

            S=conj(V1).*U1+conj(V2).*U2+conj(V3).*U3;

            M[i,j]=sum(S[:].*W[:]);
        end
    end

    return M;
end

 function curl_basis_mass_matrix(p::Array{Int64},order::Int64)
    n=length(p);
    M=Array(Complex128,(n,n));
    for i = 1 : n
        for j = 1 : n
           M[i,j]=curl_basis_inner_product(order,p[i],p[j]);
        end
    end

    return M;
end

 function curl_basis_mass_matrix(order::Int64)
    n=ordertolen(order);
    M=spzeros(Complex128,n,n);
    for i = 1 : n
        for j = 1 : n
           M[i,j]=curl_basis_inner_product(order,i,j);
        end
    end
    return M;
end

 function eval_curl_soln(
    u::Array{Complex128},
    p::Array{Int64},
    order::Int64,
    x,y,z)

    outx::Complex128=0.0;
    outy::Complex128=0.0;
    outz::Complex128=0.0;
    for i = 1 : length(u)
        (ux,uy,uz)=evaluate_curl_basis(order,p[i],x,y,z);
        outx=outx+u[i]*ux;
        outy=outy+u[i]*uy;
        outz=outz+u[i]*uz;
    end
    return (outx,outy,outz);
end


function curl_basis_vandermondes(
    order::Int64,
    p::Array{Int64},
    qx::Array{Float64},
    qy::Array{Float64},
    qz::Array{Float64},
    qw::Array{Float64})

    n=length(p);
    nq=length(qw);
    V=Array(Array{Complex128},(3,));
    V[1]=Array(Complex128,(nq,n));
    V[2]=Array(Complex128,(nq,n));
    V[3]=Array(Complex128,(nq,n));

    for i = 1 : n
        for iq = 1 : nq
            (i1,i2,i3)=ind2sub(size(qw),iq);
            (ux,uy,uz)=evaluate_curl_basis(order,p[i],qx[i1],qy[i2],qz[i3]);
            V[1][iq,i]=ux;
            V[2][iq,i]=uy;
            V[3][iq,i]=uz;
        end
    end
    return V;
end

function curl_basis_gradient_vandermondes(
    order::Int64,
    p::Array{Int64},
    qx::Array{Float64},
    qy::Array{Float64},
    qz::Array{Float64},
    qw::Array{Float64})

    n=length(p);
    nq=length(qw);
    GV  = Array(Array{Complex128},(3,3));
    for i = 1 : 3
        for j = 1 : 3
            GV[i,j]=Array(Complex128,(nq,n));
        end
    end

    for i = 1 : n
        for iq = 1 : nq
            (i1,i2,i3)=ind2sub(size(qw),iq);
            U=evaluate_curl_basis_gradient(order,p[i],qx[i1],qy[i2],qz[i3]);
            for di=1:3
                for dj=1:3
                    GV[di,dj][iq,i]=U[di,dj];
                end
            end
        end
    end
    return GV;
end
