%*********************************************************************
% Auxiliar code 
% L2 projection of a solution on a triangular mesh
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [proj_sol,new_geo]= Projection_meshes1(Solution,original_geo,new_geo)

load('Parameters.mat','N_micro')

proj_sol = Solution;
if isempty(new_geo)
    dx = 1/N_micro;
    x = -0.5+dx/2:dx:0.5-dx/2;
    [ym,xm] = meshgrid(x,x);
    
    xx = xm-dx/2;
    xx = [xx xx(:,end)];
    xx = [xx;0.5*ones(1,size(xx,2))];
    yy = ym-dx/2;
    yy = [yy;yy(end,:)];
    yy = [yy 0.5*ones(size(yy,1),1)];
    
    new_geo = struct();
    new_geo.coordinate = [xx(:),yy(:)];
    new_geo.element = delaunay(xx,yy);
    new_geo.TR = triangulation(new_geo.element,new_geo.coordinate); 
    [new_geo] =  myedge(new_geo);
end

proj_sol.Pres     = zeros(size(new_geo.element,1),1);
proj_sol.Vel      = zeros(3*new_geo.nElement,2);

elem              = pointLocation(original_geo.TR,new_geo.bari);
proj_sol.Pres(:,1)= Solution.Pres(elem);
vel_elem          = [3*(elem-1)+1; 3*(elem-1)+2; 3*(elem-1)+3];
proj_sol.Vel(:,1:2) = Solution.Vel(vel_elem,:);





