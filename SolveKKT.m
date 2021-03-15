function [x1,x2,x3]=SolveKKT(Lh,Ls,invHAT,b1,b2,b3,solution)
global A G
s=solution.s;
z=solution.z;
b1=b1-A'*b2;
g=G'*(z./s.*b3)-b1;
invHg=Lh'\(Lh\g);
x2=Ls'\(Ls\(A*invHg-b2));
x1=invHg-invHAT*x2;
x3=z./s.*(G*x1-b3);
end