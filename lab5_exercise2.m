clc;
clear;
close all;
lab5_ex2(100,10)
function [y,s] = lab5_ex2(N,M)

h = 2/N;

tim = -1:h:1;

s0 = (exp(1) - exp(-1))/(exp(1) + exp(-1))^2;

s = zeros(1,M+1);
s(1) = s0;
y = zeros(4,N+1);
y(:,1) = [(exp(1) + exp(-1))^-1 ; s0;0;1];
fun = @(t,p) [p(2);
                       -p(1)+2*(p(2)^2)/p(1) ;
                       p(4);
                       (-1 -2*(p(2)/p(1))^2)*p(3) + (4*p(2)/p(1))*p(4)];
for j = 1:M
    y(:,1) =  [(exp (1)+ exp(-1))^-1 ; s(j);0;1];
    for i = 1:N
        
        p1 = y(:,i);
        k1 = fun(tim(i),p1);
        p2 = y(:,i) + h*0.5*k1;
        k2 = fun(tim(i)+0.5*h,p2);
        p3 = y(:,i) + h*(0.5*k2);
        k3 = fun(tim(i)+0.5*h,p3);
        p4 = y(:,i) + h*k3;
        k4 = fun(tim(i)+h,p4);
        y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
    end
    s(j+1) = s(j) - (y(1,N+1)- 1/(exp(1) + exp(-1)))/y(3,N+1);

end

%disp(s);%
figure;
plot(y(1,:));
hold on
plot(y(2,:));



end