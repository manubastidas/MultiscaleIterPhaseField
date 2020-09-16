%*********************************************************************
% Initial Macro-scale solution 
%*********************************************************************
%
%***------------------------------------
%***Inputs: Parameters phase field
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Sol] = InitialSolution(Macrogeo)

global u_eq p_init
% Sol = struct();
% Sol.u = u_eq*ones(Macrogeo.nElement,1);
% Sol.p = p_init*ones(Macrogeo.nElement,1);
% 
% Sol.react(:,1) = (Sol.u/0.5).^2-1;
% 
% Sol.uCont = u_eq*ones(Macrogeo.nnodes,1);
% Sol.pCont = p_init*ones(Macrogeo.nnodes,1);
% Sol.gradPcont = zeros(Macrogeo.nnodes,2);

Sol = struct();
Sol.u = u_eq*ones(Macrogeo.nElement,1);
Sol.p = p_init*ones(Macrogeo.nElement,1);

Sol.react(:,1) = (Sol.u/0.5).^2-1;

Sol.uCont = u_eq*ones(Macrogeo.nnodes,1);
Sol.pCont = p_init*ones(Macrogeo.nnodes,1);
Sol.gradPcont = zeros(Macrogeo.nnodes,2);




