clc;
clear;
close all;
bookques2(200,1,1,5,95,5,0);
function[] = bookques2(N,T,c,d,y10,y20,y30)
%Euler's method
yE = zeros(3,N+1);
yE(:,1) = [y10;y20;y30];
h = T/N;
t = 0:h:T;

for n = 1:N
    yE(:,n+1) = yE(:,n) + h*[-c*yE(1,n)*yE(2,n);c*yE(1,n)*yE(2,n) - d*yE(2,n);d*yE(2,n)];
end  
subplot(3,1,1);
plot(t,yE,'r*--');
title('Eulers  method');
xlabel('time');
ylabel('yE');
grid on;

% Backward Euler's method
ybe = zeros(3,N+1);
ybe(:,1) = [y10;y20;y30];
h = T/N;

for n = 1:N
    y_next = ybe(:,n);
    max_iter = 50;
    tol = 1e-5;
    for k = 1:max_iter
        J = [1-h*c*y_next(2),h*c*y_next(1),0;-h*c*y_next(2),1-h*c*y_next(1)+h*d,0;0,-h*d,1];
        f = [y_next(1)-ybe(1,n)+h*c*y_next(1)*y_next(2);y_next(2)-ybe(2,n)-h*c*y_next(1)*y_next(2)+h*d*y_next(2);y_next(3)-ybe(3,n)-h*d*y_next(2)];
        y = J\f;
        y_next = y_next -y;
        diff = norm(y_next);
        if(diff<tol)
            break;
        end    

    end
    ybe(:,n+1) = y_next;

end

subplot(3,1,2);
plot(t,ybe,'b+--');
title('Backward Eulers  method');
xlabel('time');
ylabel('ybe');
grid on;

%trapezoidal method

yt = zeros(3,N+1);
yt(:,1) = [y10;y20;y30];

for n = 1:N
    
    ytn = yt(:,n);
    max_iter = 50;
    tol = 1e-5;
    for k = 1:max_iter
        J = [2+h*c*ytn(2),h*c*ytn(1),0;-h*c*ytn(2),2-h*c*ytn(1),0;0,-h*d,2];
        f = [2*ytn(1)+h*c*(yt(1,n)*yt(2,n)+ ytn(1)*ytn(2))-2*yt(1,n);2*ytn(2)-h*c*(yt(1,n)*yt(2,n) + ytn(1)*ytn(2))-2*yt(2,n);2*ytn(3)-2*yt(2,n)-h*d*(ytn(2)+yt(2,n))];
        y = J\f;
        ytn = ytn -y;
        diff = norm(ytn);
        if(diff<tol)
            break;
        end    

    end
    yt(:,n+1) = ytn;


end    
subplot(3,1,3);
plot(t,yt,'g+--');
title('Trap  method');
xlabel('time');
ylabel('yt');
grid on;

end