clc;
clear;
close all;

ques(500,20);
function [u] = ques(N,m)
h = 1/(N+1);
t = 0:h:1;

A = zeros(N,N);
for i = 1:N
    A(i,i) = -2;
end
for i = 1:N-1
    A(i,i+1) = 1;
end
for i = 2:N
    A(i,i-1) = 1;
end

u = zeros(N,1);
for i = 1:N
    u(i,1) = (i-1)/N;
end
g = zeros(N,1);
g(N,1) = -1/h^2;
f = @(u,t) 10*u^3 + 3*u + t^2;

for k = 1:m
    F = zeros(N,1);
    for idx = 1:N
        F(idx,1) = f(u(idx),idx/(N+1));
    end
    Fd = zeros(N,N);
    for idx = 1:N
        Fd(idx,idx) = 30*(u(idx)^2) + 3;
    end 
    hoo = A*u/(h*h) - F - g;
    hee = A/(h*h) - Fd;
    yoo = hee\hoo;
    u = u - yoo;
end

tim = zeros(N,1);
for idx = 1:N
    tim(idx) = idx/(N+1);
end

plot(tim,u);



end