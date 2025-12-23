clc;
clear;
N = 128;
xv = linspace(-1,3,N+1);
yv = linspace(-1,1,N+1);
h = 1.0/N;

nint = (N-1)^2;
g = @(x,y) log((x-3).^2 + (y-2).^2);


A = sparse(nint,nint);
b = zeros(nint,1);

for j = 1:N-1
    for i = 1:N-1
        k = (j-1)*(N-1) + i;
        A(k,k)= -4/h^2;
        kright = (j-1)*(N-1) + i+1;
         x_ij = xv(i+1);
         y_ij = yv(j+1);
        if(i+1 <= N-1)
            A(k,kright) = 1/h^2;
        else
            
            b(k) = b(k) - g(xv(end),y_ij)/(h*h);
        end 

        kleft = (j-1)*(N-1) + i-1;
        
        if(i-1 >= 1)
            A(k,kleft) = 1/h^2;
        else
            b(k) = b(k) - g(xv(1),y_ij)/h^2;
        end 

        kup = (j)*(N-1) + i;
        
        if(j+1 <= N-1)
            A(k,kup) = 1/h^2;
        else
            b(k) = b(k) - g(x_ij,yv(end))/h^2;
        end
        
        kdown = (j-2)*(N-1) + i;
        
        if(j-1 >=1)
            A(k,kdown) = 1/h^2;
        else
            b(k) = b(k) - g(x_ij,yv(1))/h^2;
        end

    end 
end

uvec = A\b;
% Map solution to full grid (include boundary)
U_mat = zeros(N+1, N+1);
for j = 1:(N-1)
    for i = 1:(N-1)
        k = (j-1)*(N-1) + i;
        U_mat(j+1,i+1) = uvec(k);
    end
end

% Set boundary values using g(x,y)
U_mat(1,:)   = g(xv, yv(1));       % bottom (y=0)
U_mat(end,:) = g(xv, yv(end));       % top (y=1)
U_mat(:,1)   = g(xv(1), yv);         % left (x=0)
U_mat(:,end) = g(xv(end), yv);         % right (x=1)

[X, Y] = meshgrid(xv, yv);
figure; surf(X, Y, U_mat); shading interp; colorbar;
xlabel('x'); ylabel('y'); zlabel('u(x,y)');
title('Solution with u = log((x-3)^2+(y-2)^2) on the boundary');

