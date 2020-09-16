%*********************************************************************
% Preprocess of the solution of the cell problems 
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
function [geometry] = MFEM_preProcess_Micro(geometry)

%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
% global N_macro 
% Time gamma lambda 

% Assemble matrices B and C
geometry.B = sparse(geometry.noedges, geometry.noedges);
geometry.C = sparse(geometry.noedges,geometry.nElement);
geometry.D = sparse(geometry.nElement,geometry.nElement);

for j = 1:geometry.nElement
    coord = geometry.coordinate(geometry.element(j,:),:)';
    I = full(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),geometry.element(j,[3 1 2]))));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    n = coord(:,[3,1,2])-coord(:,[2,3,1]);
    geometry.C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    geometry.D(j,j) = geometry.area(j);
    
    geometry.B(I,I)= geometry.B(I,I)+ diag(signum)*...
        stimaB(coord,eye(2))*diag(signum);
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