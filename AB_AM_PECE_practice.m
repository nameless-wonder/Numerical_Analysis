clc;
clear;
close all;
heun(0,0.005,500);
function [t, y_exact, y_pece] = heun(y0, h, N)

fun = @(t,y) -100*y + 100*t + 101;
tru = @(t) t + 1 + (y0-1)*exp(-100*t);

t = 0:h:N*h;
y_exact = tru(t);

y_pece = zeros(1,N+1);
y_pece(1) = y0;

y_pece(2) = y_pece(1) + h*fun(t(1),y_pece(1));

for n = 2:N
    pn = y_pece(n) + 0.5*h*(3*fun(t(n),y_pece(n)) - fun(t(n-1),y_pece(n-1)));
    y_pece(n+1) = y_pece(n) + h*(5*fun(t(n+1),pn) + 8*fun(t(n),y_pece(n)) - fun(t(n-1),y_pece(n-1)))/12;

end  

figure;
plot(t,y_exact,"r+--");
hold on;
plot(t,y_pece,"gO--");
grid on;






end



