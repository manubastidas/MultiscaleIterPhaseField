%*********************************************************************
% Refinement of the micro-scale mesh on the transition zone
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [geometry1,geo_full,Sol_Phasefield]= AposterioriRefinement1(Sol_Phasefieldtest,Sol_Phasefield,geometry)

load('Parameters.mat','Theta0','Theta1','lambda')

[~,coarse_geo]= Projection_meshes1(Sol_Phasefieldtest,geometry,[]);

large = Sol_Phasefieldtest.Pres >= Theta1*lambda;
small = Sol_Phasefieldtest.Pres <= 1-Theta1*lambda;
eta_T = geometry.bari(large&small,:);

geometry1  = coarse_geo;
resolution = inf;

while resolution>Theta0
    
    %% Refine mesh
    % REFINAR  = find(large&small);
    REFINAR = pointLocation(geometry1.TR,eta_T);
    REFINAR = REFINAR(~isnan(REFINAR));
        
    New_elem =  geometry1.element(REFINAR,:);
    New_points = (geometry1.coordinate(New_elem(:,[1 2 3]),:)+...
        geometry1.coordinate(New_elem(:,[2 3 1]),:))./2;
    New_points = unique(New_points,'rows');
    
%     [~,distance] = nearestNeighbor(geometry1.TR,New_points);
%     New_points = New_points(distance>Theta1,:);
%     if size(New_points,1)>1
%         pairs     = nchoosek(1:size(New_points,1),2);
%         distanceI = sqrt(sum((New_points(pairs(:,1),:)-New_points(pairs(:,2),:)).^2,2));
%         pairs     = pairs(distanceI<Theta1,1);
%         New_points(pairs,:)=[];
%     end
    
    %% To manage the periodicity
    limit = 0.5;
    New_points((New_points(:,1)== -limit),:) = [];
    New_points((New_points(:,1)==  limit),:) = [];
    New_points((New_points(:,2)== -limit),:) = [];
    New_points((New_points(:,2)==  limit),:) = [];
    
    %%
    Auxiliar = unique([[geometry1.coordinate zeros(size(geometry1.coordinate,1),1)];...
        New_points ones(size(New_points,1),1)],'rows');
    
    geo1            = struct();
    geo1.coordinate = Auxiliar(:,1:2);
    geo1.indic      = Auxiliar(:,3);
    geo1.element    = delaunay(geo1.coordinate(:,1),geo1.coordinate(:,2));
    geo1.TR         = triangulation(geo1.element,geo1.coordinate);
    
    [geo1]     = myedge(geo1);
    resolution = geo1.minsize;
    
    geometry1 = geo1;
end

%% Merge meshes
geo_full = struct();
Auxiliar = unique([geometry.coordinate geometry.indic;geometry1.coordinate geometry1.indic],'rows');
geo_full.coordinate = Auxiliar(:,1:2);
geo_full.indic      = Auxiliar(:,3);
geo_full.element    = delaunay(geo_full.coordinate(:,1),geo_full.coordinate(:,2));
geo_full.TR         = triangulation(geo_full.element,geo_full.coordinate);

[geo_full] =  myedge(geo_full);

[Sol_Phasefield]= Projection_meshes1(Sol_Phasefield,geometry,geo_full);
