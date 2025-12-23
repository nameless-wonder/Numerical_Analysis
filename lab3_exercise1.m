clc;
clear;
close all;
lab3_ex1(0 , 0.001 , 50);
function[y,yex,yem] = lab3_ex1(y0,h,N)

t = 0:h:N*h;

%tru = @(t) t + 0.01 + (y0 - 0.01)*exp(-100*t);

tru = @(t) t + 1 + (y0 - 1) * exp(-100 * t);


y = zeros(N+1,1);

yex = zeros(N+1,1);
yem = zeros(N+1,1);

fun = @(y,t) -100*y + 100*t + 101;

yem(1) = y0;
yem(2) = yem(1) + h * fun(yem(1), t(1));


yex(1) = y0;
yex(2) = tru(h);

for n = 2:N
    yem(n+1) = (yem(n) + 500*h*t(n+1)/12 + 505*h/12 + 2*h*fun(yem(n),t(n))/3 - h*fun(yem(n-1),t(n-1))/12)/(1 + 500*h/12);
    yex(n+1) = (yex(n) + 500*h*t(n+1)/12 + 505*h/12 + 2*h*fun(yex(n),t(n))/3 - h*fun(yex(n-1),t(n-1))/12)/(1 + 500*h/12);

end
y = tru(t);
    


 figure;
plot(t, y, 'ro--', 'DisplayName', 'Exact');
hold on;
plot(t, yex, 'b+--', 'DisplayName', 'True Start');
hold on;
plot(t, yem, 'g-', 'DisplayName', 'Euler Start');
xlabel('Time t');
ylabel('y(t)');
title('Comparison of Exact and Adams-Moulton Approximations');
legend;
grid on;




end
