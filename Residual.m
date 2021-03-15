function Res=Residual(solution)
global c A b G h
x=solution.x;
s=solution.s;
y=solution.y;
z=solution.z;
k=solution.k;
t=solution.t;
rx=-A'*y-G'*z-c*t;
ry=A*x-b*t;
rz=G*x+s-h*t;
rt=b'*y+h'*z+c'*x+k;
rsz=s.*z;
rtk=t*k;
Res.rx=rx;
Res.ry=ry;
Res.rz=rz;
Res.rt=rt;
Res.rsz=rsz;
Res.rtk=rtk;
end