  using Devectorize
include("incompressible_basis.jl");

  function ns_gradgrad(
    u::Array{Complex128},
    p::Array{Int64},
    order::Int64,
    quadx::Array{Float64},
    quady::Array{Float64},
    quadz::Array{Float64},
    quadw::Array{Float64})

    n=length(u);
    gradgrad=Array(Complex128,(n,));


      function outergrad(i::Int64,u,p,order,quadx,quady,quadz,quadw)
        dim=3;
        n=length(u);
        nq=length(quadw);
        gradgradq=Array(Complex128,(nq,));
       #Loop through quadrature nodes.
        for iq=1:nq
            (i1,i2,i3)=ind2sub(size(quadw),iq);
            gradv=evaluate_curl_basis_gradient(order,p[i],quadx[i1],quady[i2],quadz[i3]);
            gradu=zeros(Complex128,(dim,dim));
            for j=1:n
                tgradu=evaluate_curl_basis_gradient(order,p[j],quadx[i1],quady[i2],quadz[i3]);
                @devec gradu=gradu+u[j].*tgradu;
            end
            grad_dot_grad::Complex128 = 0.0;
            for d1 = 1 : dim
                for d2 = 1 : dim
                    grad_dot_grad=grad_dot_grad+gradu[d1,d2]*conj(gradv[d1,d2]);
                end
            end
            gradgradq[iq]=grad_dot_grad;
        end
        return sum(gradgradq.*quadw[:]); 
    end

   return Array{Complex128}(pmap(i->outergrad(i,u,p,order,quadx,quady,quadz,quadw),1:n));
end


  function ns_convection(
    u::Array{Complex128},
    p::Array{Int64},
    order::Int64,
    quadx::Array{Float64},
    quady::Array{Float64},
    quadz::Array{Float64},
    quadw::Array{Float64})


    dim=3;
    n=length(u);
    nq=length(quadw);

    conv=Array(Complex128,(n,));

    function outerconv(i::Int64,u,p,order,quadx,quady,quadz,quadw)
        n=length(u);
        nq=length(quadw);
        dim=3;
        convq=Array(Complex128,(nq,));
        #Loop through "test" degrees of freedom.
            #Loop through quadrature nodes.
            for iq=1:nq
                #Compute test   function.
                (i1,i2,i3)=ind2sub(size(quadw),iq);
                (vx,vy,vz) =evaluate_curl_basis(order,p[i],quadx[i1],quady[i2],quadz[i3]);
                gradu=zeros(Complex128,(dim,dim));
                ux=0.0*im;
                uy=0.0*im;
                uz=0.0*im;
                for j =1:n
                    #Compute gradient of solution.
                    tgradu=evaluate_curl_basis_gradient(order,p[j],quadx[i1],quady[i2],quadz[i3]);
                    @devec gradu=gradu+u[j].*tgradu;
                    (tux,tuy,tuz) =evaluate_curl_basis(order,p[j],quadx[i1],quady[i2],quadz[i3]);
                    ux=ux+u[j]*tux;
                    uy=uy+u[j]*tuy;
                    uz=uz+u[j]*tuz;
                end
                udivu=ux*gradu[:,1]+uy*gradu[:,2]+uz*gradu[:,3];
                convq[iq]=udivu[1]*conj(vx)+udivu[2]*conj(vy)+udivu[3]*conj(vz);
            end
            return sum(convq.*quadw[:]);
        end
    return Array{Complex128}(pmap(i->outerconv(i,u,p,order,quadx,quady,quadz,quadw),1:n));
end



  function ns_eval(
    RTR, #Cholesky factor of mass matrix
    u::Array{Complex128},
    p::Array{Int64},
    order::Int64,
    quadx::Array{Float64},
    quady::Array{Float64},
    quadz::Array{Float64},
    quadw::Array{Float64})

    conv=ns_convection(u,p,order,quadx,quady,quadz,quadw);
    gradgrad=ns_gradgrad(u,p,order,quadx,quady,quadz,quadw);
    return RTR\(gradgrad-conv);
end

function ns_eval(RTR,u::Array{Complex128},V,GV,W)

    dim=3;
    n=length(u);
    gg=zeros(u);
    #Compute grad-grad term.
    for i=1:dim
	for j=1:dim
	    gg[:]=gg[:]+(GV[i,j]')*(W*(GV[i,j]*u));
	end
    end

    #Compute convection term.
    co=zeros(u);
    for i=1:dim
	for j = 1 : dim
	    co[:]=co[:]+V[j]'*(W*((V[i]*u).*(GV[j,i]*u)))
	end
    end


    return RTR\(-gg-co);
end
