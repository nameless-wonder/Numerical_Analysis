clc;
clear;
close all;
lab6_exercise2(500,50);
function [u] = lab6_exercise2(N,m)

a = -1;
b=1;
h = (b-a)/(N+1);
alpha = (exp(1)+exp(-1))^-1;
beta = (exp(1)+exp(-1))^-1;
g = zeros(N,1);
g(1,1) = -alpha/(h*h);
g(N,1) = -beta/(h*h);

A = zeros(N,N);
for idx = 1:N
    A(idx,idx) = -2;
end
for idx = 2:N
    A(idx,idx-1) = 1;
end
for idx = 1:N-1
    A(idx,idx+1) = 1;
end
u = zeros(N,m+1);
for idx = 1:N
    u(idx,1) = alpha;
end
f = @(y1,y2) -y1 + 2*(y2)^2/y1;
for k = 1:m
    Fd = zeros(N,N);
    Fd(1,1) = -1 -0.5*(u(2,k) - alpha)*(u(2,k) - alpha)/(h*h*u(1,k)*u(1,k));
    for idx = 2:N-1
        Fd(idx,idx) = -1 - 0.5*((u(idx+1,k)-u(idx-1,k))/(h*u(idx,k)))^2;
    end
    Fd(N,N) = -1 -0.5*(beta - u(N-1,k))*(beta - u(N-1,k))/(h*h*u(N,k)*u(N,k));
    Fd(N,N-1) = -(beta-u(N-1,k))/(h*h*u(N,k));
    for idx = 2:N-1
        Fd(idx,idx-1) = -(u(idx+1,k)-u(idx-1,k))/(h*h*u(idx,k));
    end
    Fd(1,2) =  (u(1+1,k)-alpha)/(h*h*u(1,k));
    for idx = 2:N-1
        Fd(idx,idx+1) = (u(idx+1,k)-u(idx-1,k))/(h*h*u(idx,k));
    end
    F = zeros(N,1);
    for idx = 2:N-1
        F(idx,1) = f(u(idx,k),(u(idx+1,k) - u(idx-1,k))/(2*h));
    end
    F(1,1) = f(u(1,k),(u(2,k) - alpha)/(2*h));
    F(N,1) = f(u(N,k),(beta - u(N-1,k))/(2*h));
    J = (A/(h*h) - Fd);
    roo = (A*u(:,k)/(h*h) - F - g);
    yoo = J\roo;
    u(:,k+1) = u(:,k) - yoo;

end    
plot(u(:,m+1));
disp(u(:,m+1));


end

