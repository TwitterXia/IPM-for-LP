%% problem scale
m=100;
n=50;
p=20;
%% gengerate a problem
rng(5)
G=randn(m,n);
xfeas=randn(n,1);
h=G*xfeas+rand(m,1);
A=randn(p,n);
b=A*xfeas+rand(p,1);
c=randn(n,1);
%% solve LP
cvx_begin
cvx_solver ecos
variable x(n)
minimize(c'*x)
subject to
A*x==b;
G*x<=h;
cvx_end
%% save problem parameters
LP.m=m;
LP.n=n;
LP.p=p;
LP.c=c;
LP.A=A;
LP.b=b;
LP.G=G;
LP.h=h;
LP.xopt=x;
LP.optval=cvx_optval;
save('LP.mat','LP');