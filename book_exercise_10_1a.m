clc;
clear;
close all;
ques10(500,100)
function [Y] = ques10(N,m)

a = 0;
b = 1;
h = 1/N;
tim = 0:h:1;
y = zeros(4,N+1);
f = @(t,y1,y2,y3,y4) [y2;10*y1^3+3*y1+t^2;y4;(30*y1^2 + 3)*y3];
s = -1;
y(:,1) = [0;s;0;1];

for k = 1:m
    y(:,1) = [0;s;0;1];
    for n = 1:N
        p1 = y(:,n);
        k1 = f(tim(n),y(1,n),y(2,n),y(3,n),y(4,n));
        
        p2 = y(:,n) + h*0.5*k1;
        k2 = f(tim(n)+0.5*h,p2(1),p2(2),p2(3),p2(4));
        p3 = y(:,n) + h*(0.5*k2);
        k3 = f(tim(n)+0.5*h,p3(1),p3(2),p3(3),p3(4));
        p4 = y(:,n) + h*k3;
        k4 = f(tim(n)+h,p4(1),p4(2),p4(3),p4(4));
        y(:,n+1) = y(:,n) + h*(k1+2*k2+2*k3+k4)/6;
    end
    s = s - (y(1,n+1)-1)/(y(3,n+1));
    disp(s);
end    





end