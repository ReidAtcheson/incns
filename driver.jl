include("incompressible_navier_stokes.jl");
include("tocsv.jl");

timesteps=1000;
quadlen=20;
dt=0.0001;
order=8;
p=filter_dependent_modes(order);
len=length(p);
u=rand(Complex128,(len,));
M=curl_basis_mass_matrix(p::Array{Int64},order::Int64);
RTR=cholfact(M);

(x,y,z,w)=trapezoid_tensor_product(-pi,pi,quadlen);

vx=linspace(-pi,pi,30);
vy=linspace(-pi,pi,30);
vz=linspace(-pi,pi,30);

V=curl_basis_vandermondes(order,p,x,y,z,w);
GV=curl_basis_gradient_vandermondes(order,p,x,y,z,w);
W=Diagonal(w[:]);

fn=tu->ns_eval(RTR,tu,V,GV,W);
#RK4
for n=1:timesteps
    t=time();
    print("writing...\n");
    tocsv(vx,vy,vz,(tx,ty,tz)->eval_curl_soln(u,p,order,tx,ty,tz),"data/ns.csv.$(n)");
    print("done writing...\n");
    print("time taken to write = $(time()-t)\n");
    t=time();
    k1=fn(u);
    k2=fn(u+0.5*dt*k1);
    k3=fn(u+0.5*dt*k2);
    k4=fn(u+dt*k3);
    u=u+(dt/6.0)*(k1+2.0*k2+2.0*k3+k4);
    print("iteration $(n) done\n");
    print("iteration time = $(time()-t)\n");
end
