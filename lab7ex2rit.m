% clc;
% clear;
% 
% function[y] = lab7_2(N)
% 
% h = 2/N;
% y = zeros(N,1);
% 
% funvec = zeros(1,N);
% for n = 1:N
%     funvec{n} = @(t) t^(n-1);
% 
% end    
% 
% 
% 
% 
% 
% 
% end


clc;
clearvars;

% Problem setup
a = -1;
b = 1;
alpha = 1/(exp(1) + exp(-1));
beta = alpha;

% Functions
f = @(t,y,y1) -y + (2*(y1.^2))./y;
dfdy = @(t,y,y1) -1 - (2*(y1.^2))./(y.^2);
dfddy = @(t,y,y1) (4*y1)./y;

M = 10;
error = zeros(1,M);

% Use po = 3:M to compute errors for multiple N, or 4:4 for just N=16
for po = 7:7
    N = 2^po;
    fprintf("Current N: %d\n", N);
    h = (b-a)/N;
    t = a + h*(0:N)'; % Column vector

    % Precompute matrices
    V = t.^(0:N); % Vandermonde: V(i,j) = t(i)^(j-1)
    D = zeros(N+1, N+1);
    for j = 2:N+1
        D(:,j) = (j-1) * t.^(j-2);
    end

    % Matrix A
    A = zeros(N+1, N+1);
    A(1,:) = a.^(0:N); % u(-1)
    A(N+1,:) = b.^(0:N); % u(1)
    for i = 2:N
        for j = 3:N+1
            A(i,j) = (j-1)*(j-2)*(t(i)^(j-3));
        end
    end

    % Exact solution
    yexact = (exp(t) + exp(-t)).^-1;

    % Initial guess
    y = (alpha)* ones(N+1,1); % Matches boundary conditions

    % Newton's method
    iter = 0;
    max_iter = 100;
    while iter < max_iter
        % Compute u and u' efficiently
        sum1 = V * y;
        sum2 = D * y;

        % Residual vector
        fvec = zeros(N+1,1);
        fvec(1) = alpha;
        fvec(N+1) = beta;
        fvec(2:N) = f(t(2:N), sum1(2:N), sum2(2:N));
        H = A*y - fvec;

        % Jacobian Fj
        Fj = zeros(N+1, N+1);
        dfdy_vec = dfdy(t(2:N), sum1(2:N), sum2(2:N));
        dfddy_vec = dfddy(t(2:N), sum1(2:N), sum2(2:N));
        Fj(2:N,:) = dfdy_vec .* V(2:N,:) + dfddy_vec .* D(2:N,:);
        dH = A - Fj;

        % Check condition number
        cond_number = cond(dH);
        fprintf('Iteration %d, Condition number: %e\n', iter, cond_number);

        % Update
        y1 = y - 0.3*(dH\H);
        if norm(y1 - y) < 1e-6
            y = y1;
            break;
        end
        y = y1;
        iter = iter + 1;
    end
    if iter == max_iter
        warning('Maximum iterations reached');
    end

    % Solution
    sol = V * y;
    error(po) = max(abs(yexact - sol))/max(abs(yexact));
end

% Display errors (adjust if po loop changes)
fprintf('\nN\tError\t\tRatio\n');
for m = 3:M
    if m <= po
        fprintf('%6d\t%0.6e\t%0.2f\n', 2^m, error(m)); % ifelse(m > 3, error(m-1)/error(m), NaN)
    end
end

% Plot
plot(t, yexact, 'LineWidth', 2);
hold on;
plot(t, sol, '--', 'LineWidth', 2);
legend('Exact', 'Computed');
