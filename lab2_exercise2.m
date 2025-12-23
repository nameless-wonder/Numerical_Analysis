clc;
clear;
close all;
lab2_exercise2(250,250,1,0.5,0.1,0.02,100,10);
function [] = lab2_exercise2(N, T, a1, a2, b1, b2, y10, y20)

% f1 = @(y1,y2) (y1*(a1-b1*y2));
% f2 = @(y1,y2) (y2*(-a2+b2*y1));
h = T/N;
t = 0:h:T;
% df11 = @(y1,y2) (a1-b1*y2);
% df12 = @(y1,y2) (-b1*y1);
% df21 = @(y1,y2) (y2*b2);
% df22 = @(y1,y2) (-a2+b2*y1);
yBE  = zeros(2,N+1);
yBE(:,1) = [y10;y20];
max_iter = 5;
tol = 1e-3;

for n = 1:N
    y_next = yBE(:,n);
    k=1;
    diff = 1;
    while k<max_iter & diff>tol
        J = [(1- h*a1 + h*b1*y_next(2)), h*b1*y_next(1); -h*b2*y_next(2),1+h*a2-h*b2*y_next(1)];
        f = [y_next(1) - yBE(1,n) - h*y_next(1)*(a1-b1*y_next(2));y_next(2)-yBE(2,n)-h*y_next(2)*(-a2+b2*y_next(1))];
        y = J\f;

        y_next = y_next - y;
        diff = norm(y_next);
        k = k+1;
    end    
    yBE(:,n+1) = y_next;

end

yt = zeros(2,N+1);
yt(:,1) = [y10;y20];

maxiter = 20;
toler = 1e-4;

for n = 1:N
    ytn = yt(:,n);
    for k = 1:maxiter
        f = [ytn(1)-yt(1,n)-0.5*h*(yt(1,n)*(a1-b1*yt(2,n)))-0.5*h*ytn(1)*(a1-b1*ytn(2));ytn(2)-yt(2,n)-0.5*h*yt(2,n)*(-a2+b2*yt(1,n))-0.5*h*ytn(2)*(-a2+b2*ytn(1))];
        J = [1-0.5*h*(a1-b1*ytn(2)),0.5*h*ytn(1)*b1;-0.5*h*ytn(2)*b2,1-0.5*h*(-a2+b2*ytn(1))];
        yoo = J\f;
        ytn = ytn - yoo;
        if(norm(yoo)<toler) 
            break;
        end    
    end
    yt(:,n+1) = ytn;
end
disp(yt);
plot(t,yt(1,:),'bo--');
title("y vs time");
%xlabel("time");
%ylabel("y1");
grid on;
hold on;
plot(t,yt(2,:),'r*--');
%title("y2 vs time");
%xlabel("time");
%ylabel("y2");
legend('y1','y2');


figure;
plot(t,yBE(1,:),'bo--');
title("y1 vs time");
xlabel("time");
ylabel("y1");
grid on;

figure;
plot(t,yBE(2,:),'rx:');
title("y2 vs time");
xlabel("time");
ylabel("y2");
grid on;

figure;
plot(t,yBE(1,:),'b*');
hold on;
plot(t,yBE(2,:),'rx');
title("y1 and y2 vs t");
legend('y1','y2');

end

