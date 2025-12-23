clc;
clear;
close all;
lab5_ex(0.5,0.0005,1000);
function [y,y1,y2,y3,y4] = lab5_ex(y0,h,N)

fun = @(t,y) -100*y + 100*t + 101;
tim = 0:h:N*h;


exa = @(t) (y0-1)*exp(-100*t) + t+1;
y = exa(tim);

y1 = zeros(1,N+1);
y1(1) = y0;

y2 = zeros(1,N+1);
y2(1) = y0;

y3 = zeros(1,N+1);
y3(1) = y0;

y4 = zeros(1,N+1);
y4(1) = y0;

%heun

for i =1:N
    y1(i+1) = y1(i) + h*(fun(tim(i),y1(i))+fun(tim(i+1),y1(i)+h*fun(tim(i),y1(i))))/2;

end    



%modified euler
for i = 1:N
    
    p1 = y2(i);
    k1 = fun(tim(i),p1);
    p2 = y2(i) + 0.5*h*k1;
    k2 = fun(tim(i)+0.5*h,p2);
    y2(i+1) = y2(i) + h*k2;

end    

%heun 3stage
for i = 1:N
    p1 = y3(i);
    k1 = fun(tim(i),p1);
    p2 = y3(i) + h*0.5*k1;
    k2 = fun(tim(i)+0.5*h,p2);
    p3 = y3(i) + h*(-1*k1 + 2*k2);
    k3 = fun(tim(i)+h,p3);
    y3(i+1) = y3(i) + h*(k1 + 2*k2 + k3)/6;

end

%RK simpson 4stage
for i  = 1:N
    p1 = y4(i);
    k1 = fun(tim(i),p1);
    p2 = y4(i) + h*0.5*k1;
    k2 = fun(tim(i)+0.5*h,p2);
    p3 = y4(i) + h*(0.5*k2);
    k3 = fun(tim(i)+0.5*h,p3);
    p4 = y4(i) + h*k3;
    k4 = fun(tim(i)+h,p4);
    y4(i+1) = y4(i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;


end    

figure;
plot(tim,y,"r");
hold on;
plot(tim,y1,"g");
hold on;
plot(tim,y2,"b");
hold on;
plot(tim,y3,"y");
hold on;
plot(tim,y4,"cy");
legend;


end