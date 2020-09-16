%*********************************************************************
% Imposition of PERIODIC boundary conditions for the stokes problem
% idd for BC for the cell problem associated to the permebaility
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [SOL1,SOL2] = BoundaryStokes(nEdgeTotal,nElement,...
    edgePar,MM,E1,E2,area)

%% ZERO AVERAGE
last = [zeros(2*nEdgeTotal,1); area];
MM  = [MM last; last' 0];
E1 = [E1; 0]; E2 = [E2; 0];
neq = 2*nEdgeTotal+nElement+1;

%% periodicity
Indic_periodic = sparse([1:neq edgePar(:,2)' nEdgeTotal+edgePar(:,2)'],...
    [1:neq edgePar(:,1)' nEdgeTotal+edgePar(:,1)'],...
    [ones(neq,1);ones(2*size(edgePar,1),1)]);

noeliminar = [edgePar(:,1);nEdgeTotal+edgePar(:,1)];
eliminar   = [edgePar(:,2);nEdgeTotal+edgePar(:,2)];
Indic_periodic(:,eliminar)=[];

MM = (Indic_periodic'*MM*Indic_periodic);
E1 = (Indic_periodic'*E1);
E2 = (Indic_periodic'*E2);

x1 = MM\E1; x2 = MM\E2;

% % Copy the solution
indSolut = setdiff(1:neq,eliminar);
%
SOL1 = zeros(neq,1);
SOL2 = zeros(neq,1);
%
SOL1(indSolut,1) = x1;
SOL2(indSolut,1) = x2;
%
SOL1(eliminar,1) =  SOL1(noeliminar,1);
SOL2(eliminar,1) =  SOL2(noeliminar,1);

