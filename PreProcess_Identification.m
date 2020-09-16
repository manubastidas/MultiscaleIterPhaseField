%*********************************************************************
%*                   PreProcess_Identification                       *
%*********************************************************************
%
%This function matches the edges of external boundary of the micro scale
%mesh to impose the periodicity condition.
%
%***------------------------------------
% Manuela Bastidas - 2017.

function Micro_geo = PreProcess_Identification(Micro_geo)

edge2element = [(1:size(Micro_geo.edge2element))' Micro_geo.edge2element];
Micro_geo.DirichletEdges = edge2element((edge2element(:,5)==0),1:4);

% cordGamma: [Number edge - start point - end point]
cordGamma = Micro_geo.coordinate(Micro_geo.DirichletEdges(:,2),:);
cordGamma = [Micro_geo.DirichletEdges(:,2) cordGamma];
Dirichlet = Micro_geo.DirichletEdges;

% Boundary Left - Right
eleg1  = sortrows(cordGamma(cordGamma(:,2)==min(cordGamma(:,2)),:),3);
eleg11 = sortrows(cordGamma(cordGamma(:,2)==max(cordGamma(:,2)),:),3);

% Boundary Bottom - Top
eleg2  = sortrows(cordGamma(cordGamma(:,3)==min(cordGamma(:,3)),:),2);
eleg22 = sortrows(cordGamma(cordGamma(:,3)==max(cordGamma(:,3)),:),2);

% Pairs of nodes with the same x or y coordinates
pares1 = [eleg1(:,1) eleg11(:,1)];
pares2 = [eleg2(:,1) eleg22(:,1)];

nparX = size(pares1,1);
nparY = size(pares2,1);

% Inicializing
% EDGE 1 - EDGE 2
Micro_geo.edgePar = zeros(nparX+nparY-2,4);

%%
%X axe 
for i = 1:nparX-1
    j=i+1;
    poslado1 = (Dirichlet(:,2)==pares1(i,1)&Dirichlet(:,3)==pares1(j,1)...
        |Dirichlet(:,2)==pares1(j,1)&Dirichlet(:,3)==pares1(i,1));
    poslado2 = (Dirichlet(:,2)==pares1(i,2)&Dirichlet(:,3)==pares1(j,2)...
        |Dirichlet(:,2)==pares1(j,2)&Dirichlet(:,3)==pares1(i,2));
    % Micro_geo.edgePar -> [edge1 , edge2, element1, element2]
    Micro_geo.edgePar(i,:) = [Dirichlet(poslado1,1) Dirichlet(poslado2,1)...
        Dirichlet(poslado1,4) Dirichlet(poslado2,4)];
    
end

%%
%Y axe
for i = 1:nparY-1
    j=i+1;
    poslado1 = (Dirichlet(:,2)==pares2(i,1)&Dirichlet(:,3)==pares2(j,1)...
        |Dirichlet(:,2)==pares2(j,1)&Dirichlet(:,3)==pares2(i,1));
    poslado2 = (Dirichlet(:,2)==pares2(i,2)&Dirichlet(:,3)==pares2(j,2)...
        |Dirichlet(:,2)==pares2(j,2)&Dirichlet(:,3)==pares2(i,2));
    % Micro_geo.edgePar -> [edge1 , edge2, element1, element2]
    Micro_geo.edgePar(nparX-1+i,:) = [Dirichlet(poslado1,1) Dirichlet(poslado2,1)...
        Dirichlet(poslado1,4) Dirichlet(poslado2,4)];
end



