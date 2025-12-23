clc;
clear;
close all;
J = 1091.3*10^-6;
mu = 3986004.418 * 10^8;
r0 = 6378137;
lab4_ex2(0,0,6000,0,20,0.02,890)
function[Px, Py,Pz] = lab4_ex2(lat0, lon0,v0, gam0,phi0, h, T)
lat0 = deg2rad(lat0);
lon0 = deg2rad(lon0);
gam0 = deg2rad(gam0);
phi0 = deg2rad(phi0);
J = 1091.3*10^-6;
mu = 3986004.418 * 10^8;
r0 = 6378137;
N = T/h;
tim = 0:h:T;

Px  = zeros(1,N+1);
Py = zeros(1,N+1);
Pz = zeros(1,N+1);

Px(1) = r0*cos(lat0)*cos(lon0);
Py(1) = r0*cos(lat0)*sin(lon0);
Pz(1) = r0*sin(lat0);

Vx  = zeros(1,N+1);
Vy = zeros(1,N+1);
Vz = zeros(1,N+1);

Vini = [-sin(lat0)*cos(lon0),sin(lon0),cos(lat0)*cos(lon0);
        -sin(lat0)*sin(lon0),-cos(lon0),cos(lat0)*sin(lon0);
        cos(lat0),0,sin(lat0)]*[v0*cos(gam0)*cos(phi0);v0*cos(gam0)*sin(phi0);v0*sin(gam0)];

Vx(1)= Vini(1);
Vy(1)= Vini(2);
Vz(1)= Vini(3);
U = [Px(1); Py(1); Pz(1); Vx(1); Vy(1); Vz(1)];


for i = 1:N

    p1 = U;
    
    k1 = derivU(p1);
    p2 = U + 0.5*h*k1;
    k2 = derivU(p2);
    p3 = U + 0.5*h*k2;
    k3 = derivU(p3);
    p4 = U + h*k3;
    k4 = derivU(p4);

    U = U + h*(k1/6 + k2/3 + k3/3 + k4/6);
    Px(i+1) = U(1);
    Py(i+1) = U(2);
    Pz(i+1) = U(3);
    Vx(i+1) = U(4);
    Vy(i+1) = U(5);
    Vz(i+1) = U(6);


end    

end

function dU = derivU(U)
J = 1091.3*10^-6;
mu = 3986004.418 * 10^8;
r0 = 6378137;    
Px = U(1); Py = U(2); Pz = U(3);
Vx = U(4); Vy = U(5); Vz = U(6);
r = sqrt(Px^2 + Py^2 + Pz^2);
fact_common = (r0/r)^2;
fact_J = (3/2) * J * (r0/r)^4;

factor_xy = fact_common + fact_J * (1 - 5*(Pz/r)^2);

factor_z = fact_common + fact_J * (3 - 5*(Pz/r)^2);

% Acceleration components:
ax = - mu/(r0^2) * (Px/r) * factor_xy;
ay = - mu/(r0^2) * (Py/r) * factor_xy;
az = - mu/(r0^2) * (Pz/r) * factor_z;

% Form derivative of state vector:
dU = [Vx; Vy; Vz; ax; ay; az];

end


