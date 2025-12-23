
clc; clearvars;

a = -1; b = 1;
alpha = 1/(exp(1) + exp(-1));
beta = alpha;

f = @(t,y,y1) -y + (2*(y1.^2))./y;
dfdy = @(t,y,y1) -1 - (2*(y1.^2))./(y.^2);
dfddy = @(t,y,y1) (4*y1)./y;

error = zeros(1, 8);

for po = 3:8
    N = 2^po;
    h = 2/N;
    t = linspace(a, b, N+1)';

    % Build Lagrange basis
    L = zeros(N+1, N+1);
    Lp = zeros(N+1, N+1);
    Lpp = zeros(N+1, N+1);

    for i = 1:N+1
        Li = ones(N+1,1);
        for j = 1:N+1
            if j ~= i
                Li = Li .* (t - t(j))/(t(i) - t(j));
            end
        end
        L(:,i) = Li;

        % Derivatives via finite differences
        eps = 1e-8;
        tp = t + eps;
        tm = t - eps;
        Lip = ones(N+1,1);
        Lim = ones(N+1,1);
        for j = 1:N+1
            if j ~= i
                Lip = Lip .* (tp - t(j))/(t(i) - t(j));
                Lim = Lim .* (tm - t(j))/(t(i) - t(j));
            end
        end
        Lp(:,i) = (Lip - Lim) / (2*eps);
        Lpp(:,i) = (Lip - 2*Li + Lim) / (eps^2);
    end

    % Matrix A for 2nd derivative
    A = Lpp;
    A(1,:) = L(1,:);
    A(end,:) = L(end,:);

    % Initial guess
    y = alpha + (beta - alpha) * (t - a) / (b - a);

    iter = 0; max_iter = 100;
    while iter < max_iter
        u = L * y;
        u1 = Lp * y;

        fvec = zeros(N+1,1);
        fvec(1) = alpha;
        fvec(end) = beta;
        fvec(2:end-1) = f(t(2:end-1), u(2:end-1), u1(2:end-1));

        H = A * y - fvec;

        % Jacobian
        df1 = dfdy(t(2:end-1), u(2:end-1), u1(2:end-1));
        df2 = dfddy(t(2:end-1), u(2:end-1), u1(2:end-1));
        J = zeros(N+1, N+1);
        J(2:end-1,:) = df1 .* L(2:end-1,:) + df2 .* Lp(2:end-1,:);
        dH = A - J;

        delta = dH \ H;

        % % Line search
        % candidateLambdas = linspace(0.001, 1, 100);
        % bestLambda = candidateLambdas(1);
        % bestNorm = inf;
        % 
        % for lambda_candidate = candidateLambdas
        %     y_candidate = y - lambda_candidate * delta;
        % 
        %     u_candidate = L * y_candidate;
        %     u1_candidate = Lp * y_candidate;
        % 
        %     fvec_candidate = zeros(N+1,1);
        %     fvec_candidate(1) = alpha;
        %     fvec_candidate(end) = beta;
        %     fvec_candidate(2:N) = f(t(2:N), u_candidate(2:N), u1_candidate(2:N));
        % 
        %     H_candidate = A * y_candidate - fvec_candidate;
        % 
        %     currentNorm = norm(H_candidate);
        %     if currentNorm < bestNorm
        %         bestNorm = currentNorm;
        %         bestLambda = lambda_candidate;
        %     end
        % end

        y_new = y - bestLambda * delta;

        if norm(y_new - y) < 1e-6
            y = y_new;
            break;
        end

        y = y_new;
        iter = iter + 1;
    end

    yexact = 1 ./ (exp(t) + exp(-t));
    sol = L * y;
    error(po) = max(abs(sol - yexact)) / max(abs(yexact));
end

fprintf('\nLagrange Collocation:\nN\tError\t\tRatio\n');
for m = 3:8
    if m > 3
        fprintf('%6d\t%0.6e\t%0.2f\n', 2^m, error(m), error(m-1)/error(m));
    else
        fprintf('%6d\t%0.6e\t\n', 2^m, error(m));
    end
end

