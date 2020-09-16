%*********************************************************************
% Solution of the macro-scale solute concentration 
% 
% This code is based on:
% Boffi, D., Brezzi, F., & Fortin, M. (2013).
% Mixed finite element methods and applications (Vol. 44, pp. xiv-685).
% Heidelberg: Springer.
%
%*********************************************************************
%
function [MacroSol] = MacroProb1(MacroSol,geometry,A_Efective,...
    Phi_0,Phi_1,u_past,uiter,C,D)

% Problem MACRO1 (FOR CONCENTRATION)

global ustar Time u_down 

B2  = sparse(geometry.noedges, geometry.noedges);
BB2 = sparse(geometry.noedges,geometry.nElement);

for j = 1:geometry.nElement
    
    coord = geometry.coordinate(geometry.element(j,:),:)';
    
    I = full(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),...
        geometry.element(j,[3 1 2]))));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    %% MacroProblem 1 (FOR U)
    %     elem = sprintf('elem%i',j);
    
    B2(I,I)= B2(I,I)+ diag(signum)*stimaB(coord,inv(A_Efective(:,:,j)))*diag(signum);
    
    %     Phi_1(j,1) = Phasefield1{j}.average;
    %     Phi_0(j,1) = Phasefield0{j}.average;
    
    dofp = MacroSol.fluxp(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),...
        geometry.element(j,[3 1 2]))));
    
    BB2(I,j) = BB2(I,j) + sum(repmat(dofp,1,3).*B2(I,I))';
    
end

%% Problem MACRO1 (FOR CONCENTRATION)
Lestab = 0;
A2 = [B2 , -C   ;
    Time.dt*(C'-BB2'), Phi_1'.*(D + Lestab) ];
G = D*(Phi_0'.*u_past + Phi_1'.*ustar - Phi_0'.*ustar + Lestab.*uiter);
% A2  = [A2 last; last' 0];
G = [zeros(geometry.noedges,1);G];

if ~isempty(geometry.uDir_Gamma)
    % Dirichlet conditions
    for k = 1:size(geometry.uDir_Gamma,1)
        G(geometry.nodes2edge(geometry.uDir_Gamma(k,1),geometry.uDir_Gamma(k,2)))=...
            norm(geometry.coordinate(geometry.uDir_Gamma(k,1),:)-...
            geometry.coordinate(geometry.uDir_Gamma(k,2),:))*u_down;
    end
    
    neq = geometry.noedges+geometry.nElement;
    tmp = zeros(neq,1);
    tmp(diag(geometry.nodes2edge(geometry.uNew_Gamma(:,1),geometry.uNew_Gamma(:,2))))=...
        ones(size(diag(geometry.nodes2edge(geometry.uNew_Gamma(:,1),geometry.uNew_Gamma(:,2))),1),1);
    
    FreeEdge = find(~tmp); % + nElemnt numeration
    
    impose_dir = diag(geometry.nodes2element(geometry.uDir_Gamma(:,1),geometry.uDir_Gamma(:,2)));
    
    freePOS  = (setdiff(FreeEdge,geometry.noedges+impose_dir))';

    % SOURCE TERM - Impossing Pressure
    x = zeros(neq,1);
    x(geometry.noedges+impose_dir,1) = u_down;
    
    %     G = G-A2*x;
    % LINEAR SOLUTION (DIRECT SOLVER)
    x(freePOS,1) = A2(freePOS,freePOS)\(G(freePOS)-A2(freePOS,geometry.noedges+impose_dir)*x(geometry.noedges+impose_dir,1));
    
else
    neq = geometry.noedges+geometry.nElement;
    tmp = zeros(neq,1);
    tmp(diag(geometry.nodes2edge(geometry.Gamma(:,1),geometry.Gamma(:,2))))=...
        ones(size(diag(geometry.nodes2edge(geometry.Gamma(:,1),geometry.Gamma(:,2))),1),1);
    
    FreeEdge = find(~tmp); % + nElemnt numeration
    freePOS  = (setdiff(FreeEdge,geometry.noedges+geometry.source1))';
    
    % SOURCE TERM - Impossing Pressure
    x = zeros(neq,1);
    x(geometry.noedges+geometry.source1,1) = u_down;
%     x(geometry.noedges+geometry.source2,1) = u_up;
    
    G = G-A2*x;
    % LINEAR SOLUTION (DIRECT SOLVER)
    x(freePOS,1) = A2(freePOS,freePOS)\G(freePOS);
end

% post-process
MacroSol.u = x(geometry.noedges+1:end);
MacroSol.fluxu = x(1:geometry.noedges);
MacroSol.react = (x(geometry.noedges+1:end).^2)./(0.5^2)-1;

[MacroSol.gradU,MacroSol.gradUcont]=fluxEB(geometry.element,geometry.coordinate,...
    -MacroSol.fluxu,geometry.nodes2edge,geometry.edge2element);

%% CONTINUOS SOLUTIONS
contador = zeros(size(geometry.coordinate,1),1);
MacroSol.pCont = zeros(size(geometry.coordinate,1),1);
MacroSol.uCont = zeros(size(geometry.coordinate,1),1);

for j = 1:geometry.nElement
    MacroSol.pCont(geometry.element(j,:),1) =...
        MacroSol.pCont(geometry.element(j,:),1) + MacroSol.p(j,1);
    MacroSol.uCont(geometry.element(j,:),1) =...
        MacroSol.uCont(geometry.element(j,:),1) + MacroSol.u(j,1);
    
    contador(geometry.element(j,:),1)= contador(geometry.element(j,:),1) + ones(3,1);
end
MacroSol.pCont = MacroSol.pCont./contador;
MacroSol.uCont = MacroSol.uCont./contador;

MacroSol.Mag_gradp = sqrt(sum(MacroSol.gradPcont,2).^2);

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
