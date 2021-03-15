function solution=PredictorStep(Lh,Ls,invHAT,Res,solution)
global m c b h
rx=Res.rx;
ry=Res.ry;
rz=Res.rz;
rt=Res.rt;
rsz=Res.rsz;
rtk=Res.rtk;
s=solution.s;
z=solution.z;
k=solution.k;
t=solution.t;
%% search direction
r5=rsz;
r6=rtk;
r1=-rx;
r2=-ry;
r3=-rz+r5./z;
r4=-rt+r6/t;
[dx1,dy1,dz1]=SolveKKT(Lh,Ls,invHAT,r1,r2,r3,solution);
[dx2,dy2,dz2]=SolveKKT(Lh,Ls,invHAT,c,b,h,solution);
dt=(r4-c'*dx1-b'*dy1-h'*dz1)/(c'*dx2+b'*dy2+h'*dz2-k/t);
% dx=dx1+dt*dx2;
% dy=dy1+dt*dy2;
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
%% centering parameter
gap_p=(s+alpha*ds)'*(z+alpha*dz)+(k+alpha*dk)*(t+alpha*dt);
gap_hat=sum(rsz)+rtk;
sigma=(gap_p/gap_hat)^3;
%% corrector item
dkdt=dk*dt;
dsdz=ds.*dz;
solution.sigma=sigma;
solution.dkdt=dkdt;
solution.dsdz=dsdz;
end