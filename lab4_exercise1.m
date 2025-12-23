clc;
clear;
close all;
lab4_1(0,0.005,500);
function[y,ym1,ym2] = lab4_1(y0,h,N)

fun = @(t,y) -100*y + 100*t + 101;

ym1 = zeros(1,N+1);
ym2 = zeros(1,N+1);

ym1(1) = y0;
ym2(1) = y0;
t = 0:h:N*h;
tru = @(t) t + 1 + (y0-1)*exp(-100*t);
for n = 1:N
    ym1(n+1) = ym1(n) + h*fun(t(n),ym1(n)) + h*h*0.5*(100 -100*fun(t(n),ym1(n)));
end    

for n = 1:N
    ym2(n+1) = ym2(n) + h*fun(t(n),ym2(n)) + h*h*0.5*(100 -100*fun(t(n),ym2(n))) + h*h*h*(100*(-100) + (100*100)*fun(t(n),ym2(n)))/6;
end    

y = tru(t);

figure;
plot(t,y,"r+--");
hold on;
plot(t,ym1,"go--");
hold on ;
plot(t,ym2,"b--");
grid on;


end