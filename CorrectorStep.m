function solution=CorrectorStep(Lh,Ls,invHAT,Res,solution)
global m c A b G h qx qy qz qt
x=solution.x;
s=solution.s;
y=solution.y;
z=solution.z;
k=solution.k;
t=solution.t;
dsdz=solution.dsdz;
dkdt=solution.dkdt;
sigma=solution.sigma;
mu=solution.mu;
rx=Res.rx;
ry=Res.ry;
rz=Res.rz;
rt=Res.rt;
rsz=Res.rsz;
rtk=Res.rtk;
%% search direction
r5=rsz-sigma*mu+dsdz;
r6=rtk-sigma*mu+dkdt;
r1=-rx+sigma*mu*qx;
r2=-ry+sigma*mu*qy;
r3=-rz+sigma*mu*qz+r5./z;
r4=-rt+sigma*mu*qt+r6/t;
[dx1,dy1,dz1]=SolveKKT(Lh,Ls,invHAT,r1,r2,r3,solution);
[dx2,dy2,dz2]=SolveKKT(Lh,Ls,invHAT,c,b,h,solution);
dt=(r4-c'*dx1-b'*dy1-h'*dz1)/(c'*dx2+b'*dy2+h'*dz2-k/t);
dx=dx1+dt*dx2;
dy=dy1+dt*dy2;
dz=dz1+dt*dz2;
ds=-(r5+s.*dz)./z;
dk=-(r6+k*dt)/t;
%% step length
eta=0.989;
alpha_s=1;
alpha_z=1;
alpha_k=1;
alpha_t=1;
if any(ds<0)
    alpha_s=min(-s(ds<0)./ds(ds<0));
end
if any(dz<0)
    alpha_z=min(-z(dz<0)./dz(dz<0));
end
if dk<0
    alpha_k=-k/dk;
end
if dt<0
    alpha_t=-t/dt;
end
alpha=min([eta*alpha_s,eta*alpha_z,eta*alpha_k,eta*alpha_t,1]);
x=x+alpha*dx;
s=s+alpha*ds;
y=y+alpha*dy;
z=z+alpha*dz;
k=k+alpha*dk;
t=t+alpha*dt;
%% save solution
solution.x=x;solution.s=s;solution.y=y;solution.z=z;solution.k=k;solution.t=t;
pcost=c'*x/t;dcost=-(b'*y+h'*z)/t;
solution.pcost=pcost;solution.dcost=dcost;solution.gap=pcost-dcost;
solution.pres=norm(A*x/t-b);solution.dres=norm(A'*y/t+G'*z/t+c);
solution.mu=(s'*z+t*k)/(m+1);
solution.step=alpha;
end