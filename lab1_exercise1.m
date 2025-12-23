clc;
clear;
close all;
plot(lab1_exercise1(20));
function [E] = lab1_exercise1(J)

    A = [0,1,0;0,0,1;-2,-5,-4];
    for j = 0:J
        N = 2^j;
        h = 1/N;
    
        y = zeros(3,N+1);
        
        y(:,1) = [1;0;-1];
        for n = 1:N
            t = (n-1)*h;
            y(:,n+1) = y(:,n) + h*A*y(:,n) + h*[0;0;-4*sin(t)-2*cos(t)];
        end
        E(j+1) = abs(y(1,n+1)-cos(1));
    
    end    
    disp(E);
    for k = 1:J
        disp(E(k+1)/E(k))
    end    

end

