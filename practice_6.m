clc;
clear;
close all;
practice(500,4);
function [y,s,u] =  practice(N,m)

A = zeros(N,N);
h = 2/(N+1);
for id = 2:N-1
    A(id,id) = -2;
    A(id,id-1) = 1;
    A(id,id+1) = 1;
end
A(1,1) = -2;
A(2,2) = -2;
A(1,2) = 1;
A(2,1) = 1;
A(N,N-1) = 1;
A(N-1,N) = 1;
alpha = (exp(1) + exp(-1))^-1;
g = zeros(N,1);
g(1) = -(exp(1) + exp(-1))^-1 /(h*h);
g(N) = -(exp(1) + exp(-1))^-1 /(h*h);
u = zeros(N,1);
for lll = 1:N
    u(lll) = (exp(1) + exp(-1))^-1;
end    
f = @(y1,mor,lo) -y1 + 2*(((mor - lo)/(2*h))^2)/y1; 
f2 = @(y1,mor,lo) -1 -2*((mor - lo)/(2*h*y1))^2;
f3 = @(y1,mor,lo) (4*(mor-lo)/(2*h*y1))/(2*h);

for i = 1:m
    F = zeros(N,1);
    for id = 2:N-1
        F(id,1) = f(u(id),u(id+1),u(id-1));
    end
    F(1,1) = f(u(1),u(2),alpha);
    F(N,1) = f(u(N),alpha,u(N-1));
    Fd = zeros(N,N);
    for id  = 2:N-1
        Fd(id,id) = f2(u(id),u(id+1),u(id-1));
        Fd(id,id+1) = f3(u(id),u(id+1),u(id-1));
        Fd(id,id-1) = -f3(u(id),u(id+1),u(id-1));
    end
    Fd(1,1)= f2(u(1),u(2),alpha);
    Fd(N,N)=f2(u(N),alpha,u(N-1));
    Fd(N,N-1) = -(alpha-u(N-1))/(h*h*u(N));
    Fd(1,2) =  (u(1+1)-alpha)/(h*h*u(1));

    hoo = (A/(h*h) - Fd);
    hee = ((A*u)/(h*h)) - F - g;
    yoo = hoo\hee;
    u = u - yoo;

end   

disp(u);
plot(u);






end