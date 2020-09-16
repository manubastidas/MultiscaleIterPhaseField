%*********************************************************************
% Auxiliar geometry code to generate the micro-scale geometry 
% With grain or without grain. 
% This code allows us for different shapes inside of the mineral etc.
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [geo] = Make_microgeo

load('Parameters.mat','grain','c_grain','l_grain','N_micro')

if grain == 1
    
    % Micro Geometry
    geo            = struct();
    %  SOLVING MORE THAN ONE MACRO POINT
    [xx,yy] = meshgrid(-0.5:1/N_micro:0.5,-0.5:1/N_micro:0.5);

    x_del = abs(c_grain(1)-xx)<=l_grain&abs(c_grain(2)-yy)<=l_grain;
    y_del = abs(c_grain(1)-xx)<=l_grain&abs(c_grain(2)-yy)<=l_grain;
    xx(x_del) = [];
    yy(y_del) = [];
    kk1 = boundary([xx',yy']);
    Constrain1 = [kk1(1:end-1) kk1(2:end)];
    
    [xx_g,yy_g] = meshgrid(c_grain(1)-l_grain:1/N_micro:c_grain(1)+l_grain,...
        c_grain(2)-l_grain:1/N_micro:c_grain(2)+l_grain);
    xx_g = xx_g(:); yy_g = yy_g(:);
    kk2 = boundary([xx_g,yy_g]);
    
    xx2 = [xx xx_g(kk2(1:end-1))'];
    yy2 = [yy yy_g(kk2(1:end-1))'];
    new_kk = [length(xx)+1:length(xx2),length(xx)+1];
    
    Constrain2 = [new_kk(1:end-1)' new_kk(2:end)']; %[ind, [ind(2:end);ind(1)]];
    
    DT = delaunayTriangulation([xx2',yy2'],[Constrain1;Constrain2]);
    TF = isInterior(DT);
    geo.coordinate = DT.Points;
    geo.element    = DT.ConnectivityList(TF,:);
    geo.indic = geo.coordinate(:,1)*0;
    geo.TR = triangulation(geo.element,geo.coordinate);
    
    [geo] =  myedge(geo);
    geo.GammaInt = Constrain2;
else
    % JUST THE SQUARE
    dx = 1/N_micro;
    x = -0.5+dx/2:dx:0.5-dx/2;
    [ym,xm] = meshgrid(x,x);
    xx = xm-dx/2;
    xx = [xx xx(:,end)];
    xx = [xx;0.5*ones(1,size(xx,2))];
    yy = ym-dx/2;
    yy = [yy;yy(end,:)];
    yy = [yy 0.5*ones(size(yy,1),1)];
    
    geo = struct();
    geo.coordinate = [xx(:),yy(:)];
    geo.element = delaunay(xx,yy);
    geo.indic = xx(:)*0;
    geo.TR = triangulation(geo.element,geo.coordinate);
    
    % Microgeometry{j} features
    [geo] =  myedge(geo);
    geo.GammaInt = [];
end

