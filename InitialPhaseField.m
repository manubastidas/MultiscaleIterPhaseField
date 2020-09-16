function [Microgeometry,Init_Phasefield,Radii,Porosity,...
    Stop_geo,Stop_Phasefield,Stop_Radii,Stop_Porosity,A_Efective,K_Efective] =...
    InitialPhaseField(Macrogeo,R,R_stop)

global N_micro
%% PHASE FIELD 00 (BEFORE INITIAL)
[~,xx,yy] = CircularPhaseFieldT(N_micro,R);

j = 1;
geo00{j} = struct();
geo00{j}.coordinate = [xx(:),yy(:)];
geo00{j}.element = delaunay(xx,yy);
geo00{j}.indic = xx(:)*0;
geo00{j}.TR = triangulation(geo00{j}.element,geo00{j}.coordinate);

% Microgeometry{j} features
[geo00{j}] =  myedge(geo00{j});

%% PHASE FIELD INITIAL
% IF REFINEMENT
[Init_Phasefield{j},Microgeometry{j}]= AposterioriRefinement0(geo00{j},R);

[Microgeometry{j}]  = MFEM_preProcess_Micro(Microgeometry{j});

Porosity(j) = Init_Phasefield{j}.average;
Radii(j)    = Init_Phasefield{j}.radii;

% A_Efective = zeros(2,2,Macrogeo.nElement);
% K_Efective = zeros(2,2,Macrogeo.nElement);

% IF U_0 AND P_0 CONSTANTS
for j=2:Macrogeo.nElement
    Microgeometry{j}   = Microgeometry{1};
    Radii(j)           = Radii(1);
    Porosity(j)        = Porosity(1);
    Init_Phasefield{j} = Init_Phasefield{1};
    %         A_Efective(:,:,j)  = A_Efective(:,:,1);
    %         K_Efective(:,:,j)  = K_Efective(:,:,1);
end

%% PHASE FIELD STOP

[Stop_Phasefield0,Stop_geo]= AposterioriRefinement0(geo00{1},R_stop);
[Stop_geo]                 = MFEM_preProcess_Micro(Stop_geo);

[Stop_Phasefield,~,Stop_geo] = ...
    LschemeSolution_phaseField(Stop_geo,Stop_Phasefield0.Pres,...
    Stop_Phasefield0.Pres,0,10^-10);

Stop_Porosity = Stop_Phasefield.average;
Stop_Radii    = Stop_Phasefield.radii;

% CELL PROBLEMS -DIFFUSION & PERMEABILITY
[SolDiff1,SolDiff2,SolPerm1,SolPerm2] = ...
    Solution_cellProblems(Stop_geo,Stop_Phasefield);

% Effective tensors -DIFFUSION & PERMEABILITY
[A_Efective, K_Efective]= EfectiveTensor(Stop_geo,Stop_Phasefield.Pres,...
    SolDiff1,SolDiff2,SolPerm1,SolPerm2);

end

