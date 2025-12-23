clear all;
%close all;

M = input(' ');
for m = 1:M
    N = 2^m; 

    h = 1/N;
    u = zeros(N-1,N-1);    % for storing the numerical solution
    ex = zeros(N-1,N-1);   % for storing the exact solution

    %%
    %exact solution function
    exasol = @(x,y) ((x^2 + y^2)^2) / 16;

    %grid points
%     xv = linspace(0,1,N+1);
%     yv = linspace(0,1,N+1);
% Sir , please use above grid for plotting purposes
    xv = h*(0:N);
    yv = h*(0:N);
    


    %intializing matrix and vector
    nint = (N-1)^2;
    A = sparse(nint,nint);
    b  = zeros(nint,1);

    %boundary function and non homogeneous functions 
    g = @(x,y) ((x^2 + y^2)^2) / 16;
    f = @(x,y) (x^2 + y^2);

    %loop for mapping u(j,i) to uvec(k) making A.uvec = b system
    for j = 1:N-1
        for i=1:N-1
            xij = xv(i+1);      %true value of coordinate
            yij = yv(j+1);

            ex(j,i) = exasol(xij,yij);  %exact value at given coordinate
            k = (j-1)*(N-1) + i;


            A(k,k) = -4/(h^2);   % diagonal element are -4
            b(k) = (xij^2 + yij^2); %non homogeneous part
            
            %right point mapping
            if(i+1 <= N-1)
                kright = (j-1)*(N-1) + i+1;
                A(k,kright) = 1/(h^2);
            else
                b(k) = b(k) - (g(xv(end),yij)/(h^2));   %in case of neighbour
            end

            %left point mapping
            if(i-1>=1)
                kleft = (j-1)*(N-1) + i-1;
                A(k,kleft) = 1/(h^2);
            else
                b(k) = b(k) - (g(xv(1),yij)/(h^2));   %in case of neighbour
            end

            %up point mapping
            if(j+1 <= N-1)
                kup = (j)*(N-1) + i;
                A(k,kup) = 1/(h^2);
            else
                b(k) = b(k) - (g(xij,yv(end))/(h^2));   %in case of neighbour
            end

            %down point mapping
            if(j-1 >= 1)
                kdown = (j-2)*(N-1) + i;
                A(k,kdown) = 1/(h^2);
            else
                b(k) = b(k) - (g(xij,yv(1))/(h^2));   %in case of neighbour
            end
        end
    end

    uvec = A\b;   %solving  
    for j = 1:N-1
        for i = 1:N-1
            k = (j-1)*(N-1) + i;  %reverse mapping
            u(j,i) = uvec(k);
        end
    end


    %% Complete the code to file u (numerical solution) and 
    %% ex (the exact solution)
    %%

    e(m) = max(max(abs(u-ex)));
end

% figure;
% plot(e);
% figure;
% plot(u(1,:));
% hold on;
% plot(ex(1,:));
% legend;

for m = 1:M
    fprintf('%.2g ', e(m));    
end
fprintf('\n');

for m = 2:M
    fprintf('%.2g ', e(m-1)/e(m));    
end
fprintf('\n');
