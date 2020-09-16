%*********************************************************************
% Preprocess matrices for computing the macro-scale solution
% 
% This code is based on:
% Bahriawati, C., & Carstensen, C. (2005). 
% Three MATLAB implementations of the lowest-order Raviart-Thomas 
% MFEM with a posteriori error control. 
% Computational Methods in Applied Mathematics, 5(4), 333-359
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [C,D,geometry] = MFEM_preProcess_Macro(geometry,N_macro)

%%
%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************

% Assemble matrices B and C
% B = sparse(geometry.noedges, geometry.noedges);
C = sparse(geometry.noedges,geometry.nElement);
D = sparse(geometry.nElement,geometry.nElement);

geometry.size = 0;
geometry.bari = zeros(geometry.nElement,2);

% Imposition of the source as dirichlet conditions
% This part can be improved
N_real = [1 1];
int = 1/N_macro;
P = [ 0 0;int 0;0 int
    N_real(1) N_real(2); N_real(1)-int N_real(2); N_real(1) N_real(2)-int];
TR = triangulation([1 2  3; 4 5 6],P);
ele1 = []; ele2=[];

for j = 1:geometry.nElement
    
    coord = geometry.coordinate(geometry.element(j,:),:)';
    I = full(diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),geometry.element(j,[3 1 2]))));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    geometry.area(j)=polyarea(coord(1,:),coord(2,:));
    
    n = coord(:,[3,1,2])-coord(:,[2,3,1]);
    C(I,j) = diag(signum)*[norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]';
    
    D(j,j) = geometry.area(j);
    
    geometry.size = max(geometry.size,max([norm(n(:,1)) norm(n(:,2)) norm(n(:,3))]));
    %     B(I,I)= B(I,I)+ diag(signum)*...
    %         stimaB(coord,eye(2))*diag(signum);
    
    geometry.bari(j,:) = sum(coord,2)/3;
    
    %% Source
    % THIS PART SHOULD BE IMPROVED! 
    stbary= cartesianToBarycentric(TR,1,geometry.bari(j,:));
    if stbary(1)>0 && stbary(2)>0
        ele1 = [ele1;j];
    end
    stbary= cartesianToBarycentric(TR,2,geometry.bari(j,:));
    if stbary(1)>0 && stbary(2)>0
        ele2 = [ele2;j];
    end
end
geometry.source1 = ele1;
geometry.source2 = ele2;