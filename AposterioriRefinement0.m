%*********************************************************************
% Refinement of the micro-scale mesh on the transition zone
% INITAL CONDITION MESH 
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Phasefield,geometry]= AposterioriRefinement0(geometry,R)
% Construction first phase-field with mesh refinement

%% FINE MESH
global fine_N lambda Theta0 Theta1 gamma_par

[phi0_fine1,xx_fine,yy_fine] = CircularPhaseFieldT(fine_N,R);
[phi_x,phi_y] = gradient(phi0_fine1);

phi0_fine = phi0_fine1(:);
xx_fine   = xx_fine(:);
yy_fine   = yy_fine(:);

%% Mesh refinement
large = phi0_fine >= Theta1*lambda;
small = phi0_fine <= 1-Theta1*lambda;
eta_T = [xx_fine(large&small) yy_fine(large&small)];

resolution = 1;
while resolution>Theta0
    
    % REFINEMENT CRITERIA
    REFINAR = pointLocation(geometry.TR,eta_T);
    REFINAR = REFINAR(~isnan(REFINAR));
    %% Refine mesh
    New_elem   =  geometry.element(REFINAR,:);
    
    New_points = (geometry.coordinate(New_elem(:,[1 2 3]),:)+...
        geometry.coordinate(New_elem(:,[2 3 1]),:))./2;
    
    New_points = unique(New_points,'rows');
    
    %% To manage the periodicity
    limit = 0.5;
    New_points((New_points(:,1)== -limit),:) = [];
    New_points((New_points(:,1)==  limit),:) = [];
    New_points((New_points(:,2)== -limit),:) = [];
    New_points((New_points(:,2)==  limit),:) = [];
    
    %%
    Auxiliar = unique([[geometry.coordinate geometry.indic]; New_points ones(size(New_points,1),1)],'rows');
    
    geo1= struct();
    geo1.coordinate = Auxiliar(:,1:2);
    geo1.indic      = Auxiliar(:,3);
    geo1.element    = delaunay(geo1.coordinate(:,1),geo1.coordinate(:,2));
    geo1.TR = triangulation(geo1.element,geo1.coordinate);
    % Microgeometry{j} features
    [geo1]     = myedge(geo1);
    resolution = geo1.minsize;
    
    geometry = geo1;
end

%% PROJECTION
Phasefield.PresCont = zeros(geometry.nnodes,1);
contador            = zeros(geometry.nnodes,1);
Phasefield.Pres     = zeros(geometry.nElement,1);
Phasefield.average  = 0;

x_ref = -0.5:1/fine_N:0.5;y_ref = -0.5:1/fine_N:0.5;
for j= 1:geometry.nElement
    xbin = discretize (geometry.bari(j,1),x_ref);
    ybin = discretize (geometry.bari(j,2),y_ref);
    
    Phasefield.Pres(j) = phi0_fine1(ybin,xbin);
    Phasefield.PresCont(geometry.element(j,:)) = Phasefield.PresCont(geometry.element(j,:))+Phasefield.Pres(j);
    Phasefield.average = Phasefield.average + Phasefield.Pres(j)*geometry.area(j);
    
    contador(geometry.element(j,:)) = contador(geometry.element(j,:))+1;
    
    pos = 3*(j-1)+[1,2,3];
    Phasefield.Vel(pos,1) = phi_x(ybin,xbin);
    Phasefield.Vel(pos,2) = phi_y(ybin,xbin);
end
Phasefield.Vel = (gamma_par*lambda.^2)*Phasefield.Vel;
Phasefield.PresCont = Phasefield.PresCont./contador;
Phasefield.radii   = sqrt((1-Phasefield.average)/pi);

