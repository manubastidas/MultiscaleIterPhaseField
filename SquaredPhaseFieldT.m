%*********************************************************************
% Initializing squared phase-fields
% TEST CASE EXTRA
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [SOL_triang,geo] = SquaredPhaseFieldT(geo,l_phase)

global  N_micro c_phase
%l_phase is lenght of the rectangle inside
%c_phase is the center (usually 0.0)

M = N_micro;

% l_phase= l_phase+lambda;

% MESH
dx = 1/N_micro;
% dy = 1/M;
x = -0.5+dx/2:dx:0.5-dx/2;
% y = -0.5+dy/2:dy:0.5-dy/2;

% [ym,xm] = meshgrid(y,x);

% xx = xm-dx/2;
% xx = [xx xx(:,end)];
% xx = [xx;0.5*ones(1,size(xx,2))];
% yy = ym-dx/2;
% yy = [yy;yy(end,:)];
% yy = [yy 0.5*ones(size(yy,1),1)];

% PHASE-FIELD
phi = zeros(N_micro,M);
y_in = x<=(c_phase(1)-l_phase(1))|x>=(c_phase(1)+l_phase(1));
x_in = x<=(c_phase(2)-l_phase(2))|x>=(c_phase(2)+l_phase(2));
phi(x_in,:) = 1;
phi(:,y_in) = 1;

phi = [phi phi(:,end)];
phi = [phi; phi(end,:)];

%%
[geo]  = MFEM_preProcess_Micro(geo);

x = -0.5:1/N_micro:0.5;y = -0.5:1/N_micro:0.5;
phi_trian = zeros(size(geo.element,1),1);
for j= 1:size(geo.element,1)
    xbin = discretize (geo.bari(j,1),x);
    ybin = discretize (geo.bari(j,2),y);
    
    phi_trian(j,1) = phi(ybin,xbin);
end
[SOL_triang,~,geo] = LschemeSolution_phaseField(geo,phi_trian,phi_trian,0,1e-10);
% [SOL_triang,~,geo] = Solution_phaseField(geo,SOL_triang.Pres,0,1e-10);

% for j= 1:size(geo.element,1)
%     xbin = discretize (geo.bari(j,1),x);
%     ybin = discretize (geo.bari(j,2),y);
%     
%     phi_sqr(ybin,xbin) = SOL_triang.Pres(j);
% end

% figure,pcolor(xm,ym,phi'),shading flat,colorbar,title('gammel phi')
