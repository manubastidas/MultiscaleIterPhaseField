
function [solution,Indic_periodic,geometry] = ...
    LschemeSolution_phaseField(geometry,prev_it,prev_n,fu,tol)

%*********************************************************************
%*                                                                   *
%*                            MFEM - MACRO SOLUTION                  *
%*                                                                   *
%*********************************************************************
% global Time lambda gamma
load('Parameters.mat','Time','lambda','gamma_par','Lstab')


geometry = PreProcess_Identification(geometry);

neq = geometry.noedges+geometry.nElement;
Indic_periodic = sparse([1:neq geometry.edgePar(:,2)'],...
    [1:neq geometry.edgePar(:,1)'],[ones(neq,1);-ones(size(geometry.edgePar,1),1)]);
% end
Indic_periodic(:,geometry.edgePar(:,2))=[];

it = 1;
residual    = zeros();
residual(1) = inf;

pprev_it = prev_it;

[fplus,fminus,positions] = source_Lscheme(pprev_it,prev_n,fu);

% fp_prev = fp_aux;
% F- : FMINUS(Function concave)
% F+ : FPLUS (Function convex)
% positions = pprev_it>=x1&pprev_it<=x2; %% POSITIVE DERIVATIVE
% fplus            = zeros(size(fprev));
% fplus(positions) = fprev(positions);
% fminus           = fprev;
% fminus(positions)= 0;

% Scheme implicit in FMINUS
L_scheme  = 1/2*max(abs(4*lambda*fu+16*gamma_par),abs(4*lambda*fu-16*gamma_par));
fpminus   = ones(size(fplus))*L_scheme;
% fpminus(positions) = 0;

alpha = (Time.dt/(lambda^2));

while residual(it)>tol && it<20
    
    A = [(lambda.^2*gamma_par)^-1*geometry.B , -geometry.C,      ;
        alpha*geometry.C', geometry.D*(1+Lstab)+diag(alpha*geometry.D*fpminus)];
    b = [zeros(geometry.noedges,1);...
        geometry.D*(prev_n+alpha*fplus+Lstab*prev_it) + alpha*geometry.D*fminus + geometry.D*(alpha*fpminus.*pprev_it)];
    
    MM   = (Indic_periodic'*A*Indic_periodic);
    b    = (Indic_periodic'*b);
    x_it = MM\b;
    
    % Copy the solution
    indSolut = setdiff(1:neq,geometry.edgePar(:,2));
    sol_tot  = zeros(neq,1);
    
    sol_tot(indSolut,1)              = x_it;
    sol_tot(geometry.edgePar(:,2),1) =  -sol_tot(geometry.edgePar(:,1),1);
    
    it = it+1;
    residual(it) = norm(sol_tot(geometry.noedges+1:end)- pprev_it);
    pprev_it     = sol_tot(geometry.noedges+1:end);
    
    [fplus,fminus,positions] = source_Lscheme(pprev_it,prev_n,fu);
    
    fpminus            = ones(size(fplus))*L_scheme;
    %     positions          = pprev_it>=x1&pprev_it<=x2;
%     fpminus(positions) = 0;
    %     fplus              = zeros(size(fprev));
    %     fminus             = fprev;
    %     fplus(positions)   = fprev(positions);
    %     fminus(positions)  = 0;
end

solution.residual = residual;
solution.Pres     = sol_tot(geometry.noedges+1:end);
solution.flux     = sol_tot(1:geometry.noedges);

solution.average = sum(solution.Pres.*geometry.area');
solution.radii   = sqrt((1-solution.average)/pi);

[solution.Vel,solution.VelCont]=fluxEB(geometry.element,geometry.coordinate,...
    -solution.flux,geometry.nodes2edge,geometry.edge2element);

