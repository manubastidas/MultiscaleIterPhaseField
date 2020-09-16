%*********************************************************************
% Imposition of PERIODIC boundary conditions for the diff problem
% idd for BC for the cell problem associated to the Diffusivity
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [SOL1,SOL2] = BoundaryDiff(nEdgeTotal,nElement,...
    edgePar,A,F1,F2,area,Indic_periodic)

%% ZERO AVERAGE
last = [zeros(nEdgeTotal,1); area];
A  = [A last; last' 0];
F1 = [F1; 0]; F2 = [F2; 0];
neq = nEdgeTotal+nElement+1;

Indic_periodic = [Indic_periodic; [zeros(1,size(Indic_periodic,2)-1),1]];
Indic_periodic = [Indic_periodic [zeros(size(Indic_periodic,1)-1,1);1]];

A = (Indic_periodic'*A*Indic_periodic);
F1 = (Indic_periodic'*F1);
F2 = (Indic_periodic'*F2);

x1 = A\F1; x2 = A\F2;

% Copy the solution
indSolut = setdiff(1:neq,edgePar(:,2));
SOL1 = zeros(neq,1);
SOL1(indSolut,1) = x1;
SOL1(edgePar(:,2),1) =  -SOL1(edgePar(:,1),1);

SOL2 = zeros(neq,1);
SOL2(indSolut,1) = x2;
SOL2(edgePar(:,2),1) =  -SOL2(edgePar(:,1),1);