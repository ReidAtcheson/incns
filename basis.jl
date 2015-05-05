include("incompressible_basis.jl");
type Basis
    order::Int64;
end

type GradBasis
    order::Int64;
end

type SolutionState
    b::Basis;
    u::Array{Complex128};
end

type Exp
end

function dim(b::Basis)
    return ordertolen(b.order);
end

function SolutionState(b::Basis)
    n=dim(b);
    u=zeros(Complex128,(n,));
    return SolutionState(b,u);
end

function GradSolutionState(s::SolutionState)
    return GradSolutionState(s.b,s.u);
end

type GradSolutionState
    b::GradBasis;
    u::Array{Complex128};
end

function call(b::Basis,n::Int64,x,y,z)
    return evaluate_curl_basis(b.order,n,x,y,z);
end


function call(b::GradBasis,n::Int64,x,y,z)
    return evaluate_curl_basis_gradient(b.order,n,x,y,z);
end

function call(s::SolutionState,x,y,z)
    n=length(s.u);
    u=s.u;
    b=s.b;
    tu=map(i->b(i,x,y,z),1:n);
    tux=map( U -> U[1], tu);
    tuy=map( U -> U[2], tu);
    tuz=map( U -> U[3], tu);
   return (sum(tux.*u),sum(tuy.*u),sum(tuz.*u));
end


function call(s::GradSolutionState,x,y,z)
    n=length(s.u);
    u=s.u;
    b=s.b;
    return mapreduce(i->b(i,x,y,z)*u[i],(x,y)->x+y,1:n);
end
