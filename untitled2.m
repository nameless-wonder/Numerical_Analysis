clc;
clear;
close all;
heun(0,0.005,500);
function [t, y_exact, y_heun] = heun(y0, h, N)

fun = @(t,y) -100*y + 100*t + 101;
tru = @(t) t + 1 + (y0-1)*exp(-100*t);

t = 0:h:N*h;
y_exact = tru(t);

y_heun = zeros(N+1);

y_heun(1) = y0;

for n = 1:N
    
    y_heun(n+1) = y_heun(n) + 0.5*h*(fun(t(n),y_heun(n)) + fun(t(n+1),y_heun(n)+h*fun(t(n),y_heun(n))));

end  

figure;
plot(t,y_exact,"rO");
hold on;
plot(t,y_heun,"g+");





end



