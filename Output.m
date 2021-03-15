function Output(OUTPUT,iter,time)
fprintf('it       pcost       dcost       gap        pres       dres       k/t       mu       step       sigma\n');
fprintf('-------------------------------------------------------------------------------------------------\n');
for i=1:iter+1
    solution=OUTPUT{i};
    pcost=solution.pcost;dcost=solution.dcost;gap=solution.gap;
    pres=solution.pres;dres=solution.dres;
    k=solution.k;t=solution.t;mu=solution.mu;
    if i==1
        fprintf('%u   %4.3e   %4.3e   %1.0e   %1.0e   %1.0e   %1.0e   %1.0e   ---  ---\n',i-1,pcost,dcost,gap,pres,dres,k/t,mu);
    else
        step=solution.step;sigma=solution.sigma;
        fprintf('%u   %4.3e   %4.3e   %1.0e   %1.0e   %1.0e   %1.0e   %1.0e   %.4f   %1.0e\n',...
        i-1,pcost,dcost,gap,pres,dres,k/t,mu,step,sigma);
    end
end
fprintf('-------------------------------------------------------------------------------------------------\n');
if k/t<=1e-6
    fprintf('solved');
    fprintf('the optimal value: %.6f\n',pcost);
elseif pcost<0
    disp('unbounded');
elseif dcost>0
    disp('infesible');
else
    disp('failed');
end
fprintf('elpased time: %.4f\n',time);
end