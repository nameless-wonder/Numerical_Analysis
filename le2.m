clc;
clear all;
close all;
T = 1.0;

M = 12;
error = zeros(1:M);

for m = 3:M
    N = 2^m;
    h = T/N;

    y = zeros(1,N+1);   % this array holds the exact solution
    yAB = zeros(1,N+1); % ... holds the numerical solution

    %%
    % Fill the arrays y and yAB appropriately
    %3rd order taylor series required
    tim = 0:h:T;

    tru = @(t) (1/5 - 2*exp(8*t/15)/5) / (1 - 6*exp(8*t/15)/5);
    y = tru(tim);
    yAB(1) = 1;
    
    fun = @(y) (-4*y*y + 32*y/15 -4/15);


    for n = 1:4
        
        yAB(n+1) = yAB(n) + h*(-4*yAB(n)*yAB(n) + 32*yAB(n)/15 - 4/15) + h*h*0.5*((-4*yAB(n)*yAB(n) + 32*yAB(n)/15 - 4/15)*(-8*yAB(n) + 32/15)) + h*h*h*(-8*(-4*yAB(n)*yAB(n) + 32*yAB(n)/15 - 4/15) + (-8*yAB(n) + 32/15)*(-8*yAB(n) + 32/15)*(-4*yAB(n)*yAB(n) + 32*yAB(n)/15 - 4/15))/6 + h*h*h*h*(0);

    end 
    
    for n = 5:N
        
        yAB(n+1) = yAB(n) + h*(1901*fun(yAB(n)) - 2774*fun(yAB(n-1)) + 2616*fun(yAB(n-2)) -1274*fun(yAB(n-3)) + 251*(fun(yAB(n-4))))/720;

    end    


    %AB method 3 step



    
    % Do not change the template below this line
    %%


    error(m) = max(abs(y-yAB))/max(abs(y));
end

fprintf('\n');
for m = 3:M
    fprintf('%6d \t %0.6e \t %0.2f\n', ...
    2^m, error(m), error(m-1)/error(m));
end

plot([0:N]*h,y,'b++','LineWidth',2);
hold on;
plot([0:N]*h,yAB,'r--','LineWidth',2);
legend('y','yAB');