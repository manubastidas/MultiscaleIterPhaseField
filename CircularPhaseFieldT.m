%*********************************************************************
% Initial phase-field shape (circular) 
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020 (CARINA BRINGEDAL)
% Hasselt University, Belgium

function [phi,xx,yy] = CircularPhaseFieldT(N,R)
% R is wanted radius for your initial condition
% lambda is the width of the diffuse interface
% N,M are how fine grid you want it to be on (I usually use N=M, and
% N>6/lambda)

global lambda

M = N;
NN = 10^3;
dr = 1/NN;
alpha = 10;

r = (dr:dr:1)';

dt = dr;

R0 = R+lambda; % Initiall radius
f = @(r,R) 1./(1+exp(-4/lambda*(r-R)));
phi0 = f(r,R0);
phigg = f(r,R);

dx = 1/N;
dy = 1/M;

x = -0.5+dx/2:dx:0.5-dx/2;
y = -0.5+dy/2:dy:0.5-dy/2;

[ym,xm] = meshgrid(y,x);

rm = sqrt(xm.^2+ym.^2);

phig = zeros(N,M);
for i = 1:N
    for j = 1:M
        rr = rm(i,j);
        [~,J] = min(abs(rr-r));
        phig(i,j) = phigg(J);
    end
end

ex = ones(NN,1);
SS = spdiags([ex ex ex],-1:1,NN,NN);


while R0-R > 0 % While Radius too large, make timestep
%    phi0 = phi0 - dt/(alpha*xi^2)*12*phi0.*(1-phi0).*(1-2*phi0) + dt/alpha*3/4*RadialLaplace(phi0,r,dr);
    phis = phi0; % initiall guess for Newton
    F = @(t,phi) FLRphi_vektor(phi,phi0,r,dr,dt,alpha,lambda);
    J = @(phi) numjac(F,0,phi,F(0,phi0),ones(NN,1),[],0,SS,[]);
    kX = 0;
    phis = phis - J(phis)\F(0,phis);
    B = F(0,phis);
    while norm(B) > NN*10^(-14) && kX < 20
        phis = phis - J(phis)\B;
        kX = kX+1;
        B = F(0,phis);
    end
    phi0 = phis; % endeleg phi - gi til neste tidssteg
    [~,I] = min(abs(phi0-0.5));
    R0 = r(I);
end


phi = zeros(N,M);
for i = 1:N
    for j = 1:M
        rr = rm(i,j);
        [~,J] = min(abs(rr-r));
        phi(i,j) = phi0(J);
    end
end

dx = 1/N;
xx = xm-dx/2;
xx = [xx xx(:,end)];
xx = [xx;0.5*ones(1,size(xx,2))];
yy = ym-dx/2;
yy = [yy;yy(end,:)];
yy = [yy 0.5*ones(size(yy,1),1)];
%
phi = [phi phi(:,end)];
phi = [phi; phi(end,:)];

% figure,pcolor(xm,ym,phig'),shading flat,colorbar,title('gammel phi')
% figure,pcolor(xm,ym,phi'),shading flat,colorbar,title('phi')
% figure,pcolor(xm,ym,phig'-phi'),shading flat,colorbar,title('differanse')
% 
% figure,plot(r,phi0,r,f(r,R),'linewidth',2),legend('Ny phi','gammel phi')
end 

function F = FLRphi_vektor(phi,phi0,r,dr,dt,alpha,xi)
% function for Newton to use in implicit time stepping

F = alpha*xi^2*(phi-phi0)*dr^2 + 12*phi.*(1-phi).*(1-2*phi)*dr^2*dt - 3/4*xi^2*RadialLaplace(phi,r,dr)*dr^2*dt; % + 4*alpha*xi*phi.*(1-phi)*(-1)*dr^2*dt;
end
function LR = RadialLaplace(phi,r,dr)
% Radial Laplace discretization of phi

N = length(phi);

phiext = zeros(N+2,1);
phiext(1) = 0;
% phiext(1) = phi(1);
phiext(2:N+1) = phi;
% phiext(N+2) = phi(N);
phiext(N+2) = 1;

LR = 1./(r*dr).*((r+dr/2).*(phiext(3:N+2)-phiext(2:N+1))/dr - (r-dr/2).*(phiext(2:N+1)-phiext(1:N))/dr);
end

