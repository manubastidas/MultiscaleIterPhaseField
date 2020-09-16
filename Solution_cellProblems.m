%*********************************************************************
% MFEM Solution of the cell problems 
% 
% This code is based on:
% Bahriawati, C., & Carstensen, C. (2005). 
% Three MATLAB implementations of the lowest-order Raviart-Thomas 
% MFEM with a posteriori error control. 
% Computational Methods in Applied Mathematics, 5(4), 333-359
% and 
% Boffi, D., Brezzi, F., & Fortin, M. (2013).
% Mixed finite element methods and applications (Vol. 44, pp. xiv-685).
% Heidelberg: Springer.
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [SolDiff1,SolDiff2,SolPerm1,SolPerm2] = Solution_cellProblems(geometry,Sol_Phasefield)

%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
% global gamma lambda Ka n_param delta
load('Parameters.mat','gamma_par', 'lambda', 'Ka', 'n_param', 'delta')

p_prev   = Sol_Phasefield.Pres;
geometry = PreProcess_Identification(geometry);

neq            = geometry.noedges+geometry.nElement;
Indic_periodic = sparse([1:neq geometry.edgePar(:,2)'],...
    [1:neq geometry.edgePar(:,1)'],[ones(neq,1);-ones(size(geometry.edgePar,1),1)]);

Indic_periodic(:,geometry.edgePar(:,2)) = [];

%% For diffusion
alpha = (gamma_par.*lambda^2)^(-1);
B     = sparse(geometry.noedges, geometry.noedges);
grad  = zeros(geometry.nElement,2);

%% For permeability
H  = sparse(geometry.noedges, geometry.noedges);
Mx = sparse(geometry.nElement,geometry.noedges);
My = sparse(geometry.nElement,geometry.noedges);
K  = sparse(geometry.noedges, geometry.noedges);
E  = sparse(geometry.noedges,1);

p_prev = p_prev + delta;
Sol_phaseD = p_prev;

f_phi = (Ka/lambda);
f_phi = f_phi.*(((1-p_prev)*n_param)./(p_prev+n_param));
f_phi = f_phi.*(1./(Sol_phaseD).^2);

for j = 1:geometry.nElement
    coord = geometry.coordinate(geometry.element(j,:),:)';
    x     = coord(1,:); y = coord(2,:);
    
    I = full(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),...
        geometry.element(j,[3 1 2]))));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    %% Difussion Problem
    B(I,I) = B(I,I)+ (p_prev(j))^(-1)*diag(signum)*...
        stimaB(coord,eye(2))*diag(signum);
    
    % % New RHS
    grad(j,1) = sum(Sol_Phasefield.Vel(3*(j-1)+[1,2,3],1))./3;
    grad(j,2) = sum(Sol_Phasefield.Vel(3*(j-1)+[1,2,3],2))./3;
    
    %% StokesProblem
    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/geometry.area(j);
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/geometry.area(j);
    
    Sx = [-b(1)+b(2)+b(3); b(1)-b(2)+b(3); b(1)+b(2)-b(3)];
    Sy = [-c(1)+c(2)+c(3); c(1)-c(2)+c(3); c(1)+c(2)-c(3)];
    HK = (Sx*Sx'+Sy*Sy')*geometry.area(j);
    
    H(I,I)   = H(I,I)+ HK;
    Mx(j,I)  = Mx(j,I)+Sx'*geometry.area(j);
    My(j,I)  = My(j,I)+Sy'*geometry.area(j);
    K(I,I)   = K(I,I)+ f_phi(j)*diag(ones(3,1)*1/3*geometry.area(j));
    E(I)     = E(I) + (1/3)*geometry.area(j);
end

%% Diffusion
A  = [B , geometry.C   ;
    geometry.C', sparse(geometry.nElement,geometry.nElement) ];
F1 = [sparse(geometry.noedges,1); alpha.*geometry.D*grad(:,1)];
F2 = [sparse(geometry.noedges,1); alpha.*geometry.D*grad(:,2)];

[sol_tot1,sol_tot2] = BoundaryDiff(geometry.noedges,geometry.nElement,...
    geometry.edgePar,A,F1,F2,geometry.area',Indic_periodic);

% % SAVE THE SOLUTION
SolDiff1.Pres = sol_tot1(end-geometry.nElement:end-1);
SolDiff1.flux = sol_tot1(1:end-geometry.nElement-1);

SolDiff2.Pres = sol_tot2(end-geometry.nElement:end-1);
SolDiff2.flux = sol_tot2(1:end-geometry.nElement-1);

[SolDiff1.Vel,SolDiff1.VelCont] = fluxEB(geometry.element,geometry.coordinate,...
    -SolDiff1.flux,geometry.nodes2edge,geometry.edge2element);
[SolDiff2.Vel,SolDiff2.VelCont] = fluxEB(geometry.element,geometry.coordinate,...
    -SolDiff2.flux,geometry.nodes2edge,geometry.edge2element);


%% Permeability
MM = [H+K sparse(geometry.noedges,geometry.noedges) Mx';
    sparse(geometry.noedges,geometry.noedges) H+K   My';
    Mx My sparse(geometry.nElement,geometry.nElement)'];

E1 = [E;sparse(geometry.noedges,1);sparse(geometry.nElement,1)];
E2 = [sparse(geometry.noedges,1);E;sparse(geometry.nElement,1)];

% SOLUTION +BOUNDARY +ZERO AVER
[SOL1, SOL2] = BoundaryStokes(geometry.noedges,geometry.nElement,...
    geometry.edgePar,MM,E1,E2,geometry.area');

% save SOLUTION
% in thi case flux means the values over the edges
SolPerm1.flux = [SOL1(1:geometry.noedges),SOL1(1+geometry.noedges:2*geometry.noedges)];
SolPerm1.Pres = SOL1(2*geometry.noedges+1:end-1);
SolPerm2.flux = [SOL2(1:geometry.noedges),SOL2(1+geometry.noedges:2*geometry.noedges)];
SolPerm2.Pres = SOL2(2*geometry.noedges+1:end-1);

for j = 1:geometry.nElement
    %% Velocity STOKES
    nodes    = [geometry.element(j,:)' geometry.element(j,[2 3 1])'];
    edgesele = (diag(geometry.nodes2edge(nodes(:,1),nodes(:,2))));
    
    SolPerm1.Vel(j,:) = sum(SolPerm1.flux(edgesele,:))/3;
    SolPerm2.Vel(j,:) = sum(SolPerm2.flux(edgesele,:))/3;
    
    % Velocity Diffusivity
%     SolDiff1.Vel(j,:) = sum(SolDiff1.flux(edgesele,:))/3;
%     SolDiff2.Vel(j,:) = sum(SolDiff2.flux(edgesele,:))/3;
end

end 
function B=stimaB(coord,A)

%N=coord(:)*ones(1,3)-[coord;coord;coord];

N=coord(:)*ones(1,3)-repmat(coord,3,1);

D=diag([norm(N([5,6],2)) norm(N([1,2],3)) norm(N([1,2],2))]);

M=spdiags([ones(6,1),ones(6,1),2*ones(6,1),ones(6,1),ones(6,1)],...
          [-4,-2,0,2,4],6,6);

N_aux = (repmat(diag(A),3,3).*N)';    
B = D*N_aux*M*N*D/(24*det([1,1,1;coord])); 
end

