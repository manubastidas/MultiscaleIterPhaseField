%*********************************************************************
%*                                                                   *
%*              EFECTIVE Diffusion tensor A _ Upscaled               *
%*                                                                   *
%*********************************************************************
%
% This functions compute the efective difussion tensor if the original
% tensor does not dependent of the macro scale.
%
%***------------------------------------
%***Inputs: Solution of the micro cell problems.
%
%***------------------------------------
% Manuela Bastidas - 2017.

function [A_Efective, K_Efective] = EfectiveTensor(geometry,phase,...
    diffsol1,diffsol2,...
    permsol1,permsol2)

% global delta
load('Parameters.mat','delta')

Vel1= diffsol1.flux; Vel2= diffsol2.flux;

element      = geometry.element;
coordinate   = geometry.coordinate;

% AreaSum    = 0;
A_Efective = zeros(2,2);
K_Efective = zeros(2,2);
for j = 1:geometry.nElement
    % Coord (x;y) triangle vertex
    coord = coordinate(element(j,:),:)';
    
    %     coord egdes : [1:2] inicial , [3:4] final
    p      = [coord(:,[2 3 1]);coord(:,[3 1 2])];
    %     Length edge
    le     = [norm(p(1:2,1)-p(3:4,1)) norm(p(1:2,2)-p(3:4,2))...
        norm(p(1:2,3)-p(3:4,3))];
    
    I = diag(geometry.nodes2edge(geometry.element(j,[2 3 1]),...
        geometry.element(j,[3 1 2])));
    signum = ones(1,3);
    signum((j==geometry.edge2element(I,4)))= -1;
    
    N =geometry.bari(j,:)'*ones(1,3)-coord;
%     P1=N*1/2*diag(signum)*diag(le)*Vel1(I);
%     P2=N*1/2*diag(signum)*diag(le)*Vel2(I);
    
        P1 = sum(diffsol1.Vel(3*(j-1)+[1,2,3],:))/3;
        P2 = sum(diffsol2.Vel(3*(j-1)+[1,2,3],:))/3;
    A_Efective = A_Efective + (phase(j)+delta)*eye(2)*geometry.area(j) +[P1' P2']*geometry.area(j);
    
    Zeta = (phase(j)/(phase(j)+delta))*[permsol1.Vel(j,:)' permsol2.Vel(j,:)'];
    K_Efective = K_Efective + geometry.area(j)*Zeta;
    
    % A_EFECTIVE = INT_P (e_j+w^j)*e_i
    %     AreaSum = AreaSum + area;
end
% A_Efective =  A_Efective./sum(geometry.area);