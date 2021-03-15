Main
time=zeros(1,500);
for i=1:500
    tic,Main,time(i)=toc;
end
fprintf('average solve time is: %.4f\n',mean(time));