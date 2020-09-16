%*********************************************************************
% Assign nodes to the bottom line to compute the Dirichlet case
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Macrogeo,Micro_paste] = DirichletConditions(Macrogeo)

global N_macro

%% Concentration boundary (MACRO 1)
pos1 = find(Macrogeo.coordinate(Macrogeo.Gamma(:,1),1)==0);
pos2 = find(Macrogeo.coordinate(Macrogeo.Gamma(:,2),1)==0);
pos  = intersect(pos1,pos2);

Macrogeo.uDir_Gamma = Macrogeo.Gamma(pos,:);
pos_new            = setdiff(1:size(Macrogeo.Gamma,1),pos);
Macrogeo.uNew_Gamma = Macrogeo.Gamma(pos_new,:);

%% Pressure boundary (MACRO 2)
pos1 = find(Macrogeo.coordinate(Macrogeo.Gamma(:,1),1)==1);
pos2 = find(Macrogeo.coordinate(Macrogeo.Gamma(:,2),1)==1);
pos  = [pos;intersect(pos1,pos2)];

Macrogeo.pDir_Gamma = Macrogeo.Gamma(pos,:);
pos_new            = setdiff(1:size(Macrogeo.Gamma,1),pos);
Macrogeo.pNew_Gamma = Macrogeo.Gamma(pos_new,:);

%% EXAMPLE 1D
Micro_Solve = find(Macrogeo.bari(:,2)<1/N_macro);
Micro_paste = zeros(length(Micro_Solve),N_macro/2);
for j = 1:size(Micro_Solve,1)
    elem = Micro_Solve(j);
    pos_j = find(Macrogeo.bari(:,1)==Macrogeo.bari(elem,1));
    Micro_paste(j,:) = [elem setdiff(pos_j',elem)];
end
