%*********************************************************************
% Solution of the macro-scale flow problem
% 
% This code is based on:
% Bahriawati, C., & Carstensen, C. (2005). 
% Three MATLAB implementations of the lowest-order Raviart-Thomas 
% MFEM with a posteriori error control. 
% Computational Methods in Applied Mathematics, 5(4), 333-359
%
%*********************************************************************

function [MacroSol] = MacroProb2(geometry,K_Efective,C)

% Problem MACRO2 (FOR Pressure)
global p_up p_down

B1    = sparse(geometry.noedges, geometry.noedges);

for j = 1:geometry.nElement
    
    %% MacroProblem 2 (FOR P)
    coord = geometry.coordinate(geometry.element(j,:),:)';
    
    I = full(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),...
        geometry.element(j,[3 1 2]))));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    % Matrix depending on Kperm
    B1(I,I)= B1(I,I)+ diag(signum)*stimaB(coord,inv(K_Efective(:,:,j)))*diag(signum);
end

%% Problem MACRO2 (FOR Pressure)
A = [B1 , C   ;
    C', sparse(geometry.nElement,geometry.nElement) ];

neq = geometry.noedges+geometry.nElement;
F = zeros(neq,1);

if ~isempty(geometry.pDir_Gamma)
    % Dirichlet conditions
    for k = 1:size(geometry.pDir_Gamma,1)/2
        F(geometry.nodes2edge(geometry.pDir_Gamma(k,1),geometry.pDir_Gamma(k,2)))=...
            norm(geometry.coordinate(geometry.pDir_Gamma(k,1),:)-...
            geometry.coordinate(geometry.pDir_Gamma(k,2),:))*p_down;
    end
    for k = (1+size(geometry.pDir_Gamma,1)/2):size(geometry.pDir_Gamma,1)
        F(geometry.nodes2edge(geometry.pDir_Gamma(k,1),geometry.pDir_Gamma(k,2)))=...
            norm(geometry.coordinate(geometry.pDir_Gamma(k,1),:)-...
            geometry.coordinate(geometry.pDir_Gamma(k,2),:))*p_up;
    end
    
    neq = geometry.noedges+geometry.nElement;
    tmp = zeros(neq,1);
    tmp(diag(geometry.nodes2edge(geometry.pNew_Gamma(:,1),geometry.pNew_Gamma(:,2))))=...
        ones(size(diag(geometry.nodes2edge(geometry.pNew_Gamma(:,1),geometry.pNew_Gamma(:,2))),1),1);
    
    FreeEdge = find(~tmp); % + nElemnt numeration
    
    impose_dir1 = diag(geometry.nodes2element...
        (geometry.uDir_Gamma(1:size(geometry.pDir_Gamma,1)/2,1),...
        geometry.uDir_Gamma(1:size(geometry.pDir_Gamma,1)/2,2)));
    impose_dir2 = diag(geometry.nodes2element...
        (geometry.uDir_Gamma(size(geometry.pDir_Gamma,1)/2+1:end,1),...
        geometry.uDir_Gamma(size(geometry.pDir_Gamma,1)/2+1:end,2)));

    freePOS1  = (setdiff(FreeEdge,geometry.noedges+impose_dir1))';
    freePOS   = (setdiff(freePOS1,geometry.noedges+impose_dir2))';

    % SOURCE TERM - Impossing Pressure
    x = zeros(neq,1);
    x(geometry.noedges+impose_dir1,1) = p_down;
    x(geometry.noedges+impose_dir2,1) = p_up;
    
    F = F-A*x;
    % LINEAR SOLUTION (DIRECT SOLVER)
    x(freePOS,1) = A(freePOS,freePOS)\F(freePOS);
    
else
    tmp =zeros(neq,1);
    tmp(diag(geometry.nodes2edge(geometry.Gamma(:,1),geometry.Gamma(:,2))))=...
        ones(size(diag(geometry.nodes2edge(geometry.Gamma(:,1),geometry.Gamma(:,2))),1),1);
    
    FreeEdge = find(~tmp); % + nElemnt numeration
    freePOS  = (setdiff(FreeEdge,geometry.noedges+[geometry.source1;geometry.source2]))';
    
    % SOURCE TERM - Impossing Pressure
    x=zeros(neq,1);
    x(geometry.noedges+geometry.source1,1) = p_down;
    x(geometry.noedges+geometry.source2,1) = p_up;
    
    F = F-A*x;
    % LINEAR SOLUTION (DIRECT SOLVER)
    x(freePOS,1) = A(freePOS,freePOS)\F(freePOS);
end
% Post-Process
MacroSol.p     = x(geometry.noedges+1:end);
MacroSol.fluxp = x(1:geometry.noedges);

[MacroSol.gradP,MacroSol.gradPcont]=fluxEB(geometry.element,geometry.coordinate,...
    -MacroSol.fluxp,geometry.nodes2edge,geometry.edge2element);

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
