load LP.mat
%% 
tic
global m n p c A b G h qx qy qz qt
m=LP.m;
n=LP.n;
p=LP.p;
c=LP.c;
A=LP.A;
b=LP.b;
G=LP.G;
h=LP.h;
%% IMP paramerters 
Kmax=30;
iter=0;
tol=1e-10;
%% Initialize
solution=InitialPoint('trival');
Res=Residual(solution);
mu=solution.mu;
qx=Res.rx/mu;
qy=Res.ry/mu;
qz=Res.rz/mu;
qt=Res.rt/mu;
gap=mu*(m+1);
OUTPUT=cell(1,Kmax);
OUTPUT{1}=solution;
%% 
while(gap>tol)
    if iter>Kmax
        error('something must be wrong......');
    end
    iter=iter+1;
    solution.iter=iter;
    z=solution.z;s=solution.s;
    H=G'*diag(z./s)*G+A'*A;
    Lh=chol(H,'lower');
    invHAT=Lh'\(Lh\A');
    S=A*invHAT;
    Ls=chol(S,'lower');
    solution=PredictorStep(Lh,Ls,invHAT,Res,solution);
    solution=CorrectorStep(Lh,Ls,invHAT,Res,solution);
    OUTPUT{iter+1}=solution;
    Res=Residual(solution);
    gap=solution.mu*(m+1);
end
time=toc;
Output(OUTPUT,iter,time)