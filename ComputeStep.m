function alpha=ComputeStep(Solution)
s=Solution.s;
z=Solution.z;
k=Solution.k;
t=Solution.t;
ds=Solution.ds;
dz=Solution.dz;
dk=Solution.dk;
dt=Solution.dt;
%% compute the longest step one can take
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
end