clc;
clear;

alpha = 1/(exp(1) + exp(-1));
beta = alpha ;

a = -1;
b = 1;
error = zeros(1,10);
M=10;
f = @(t,y,y1) -y + (2*(y1.^2))./y;
dfdy = @(t,y,y1) -1 - 2*((y1.^2)./(y.^2));
dfddy = @(t,y,y1) 4*(y1)./(y);

for po = 1:5
    N = 2^po;
    h = 2/(N);
    t = a + h*(0:N)';

    V = t.^(0:N);
    D = zeros(N+1,N+1);

    for j = 2:N+1
        D(:,j) = (j-1)*t.^(j-2);
    
    end 

    A = zeros(N+1,N+1);

    A(1,:) = a.^(0:N);
    A(N+1,:) = b.^(0:N);

    for i = 2:N
        for j = 3:N+1
            A(i,j) = (j-2)*(j-1)*(t(i)^(j-3));
        end    
    end

    yexact = (exp(t) + exp(-t)).^-1;

    y = alpha*ones(N+1,1);

    iter =0;
    maxiter = 50;
    while iter < maxiter
        
        sum1 = V*y;
        sum2 = D*y;
        fvec = zeros(N+1,1);
        fvec(1) = alpha;
        fvec(N+1) = beta;

        fvec(2:N) = f(t(2:N),sum1(2:N),sum2(2:N));
        H = A*y - fvec;

        Fj = zeros(N+1,N+1);
        dfdy_vec = dfdy(t(2:N),sum1(2:N),sum2(2:N));
        dfddy_vec = dfddy(t(2:N),sum1(2:N),sum2(2:N));
        Fj(2:N,:) = dfdy_vec.*V(2:N,:) + dfddy_vec.*D(2:N,:);

        dH = A - Fj;

        upd = 0.5*(dH\H);

        if(norm(upd)<1e-6)
            break;
        else
            y = y-upd;
            iter = iter +1;
        end
    end 
    if iter == maxiter
        warning('Maximum iterations reached\n');
    end

    % Solution
    sol = V * y;
    error(po) = max(abs(yexact - sol))/max(abs(yexact));
    for m = 3:M
        if m <= 10
            fprintf('%6d\t%0.6e\t%0.2f\n', 2^m, error(m)); % ifelse(m > 3, error(m-1)/error(m), NaN)
        end
    end
end



% Plot
plot(t, yexact, 'LineWidth', 2);
hold on;
plot(t, sol, '--', 'LineWidth', 2);
legend('Exact', 'Computed');


  
