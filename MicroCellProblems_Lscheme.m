%*********************************************************************
% Solution of the non-linear problems and upgrading all the micro-scale
% solutions
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [Radii,Porosity,Sol_Phasefield,Sol_Phasefield_n,A_Efective,...
    K_Efective,geometry1] = MicroCellProblems_Lscheme...
    (geometry,Macrogeo,react,Sol_Phasefield_iter,Sol_Phasefield_n,...
    A_Efective0,K_Efective0,stopindic,Stop_phasefield,Stop_geometry,...
    Stop_A_Efective,Stop_K_Efective,iteracion,Micro_paste,Active,InActive,simplycopy)

%% MICRO SOLVER
global in_tol Dirichlet N_macro
% Stop_R

nElement   = Macrogeo.nElement;
A_Efective = zeros(2,2,Macrogeo.nElement);
K_Efective = zeros(2,2,Macrogeo.nElement);
Sol_Phasefield = cell(nElement,1);

geometry1 = cell(nElement,1);
Radii    = zeros(1,nElement);
Porosity = zeros(1,nElement);

stopradii_loop = Stop_phasefield.radii;
in_tol_loop    = in_tol;

if Dirichlet == 1
    To_solve = Micro_paste(:,1);
else
    To_solve = [];
end
pctRunOnAll warning off
parfor j=1:nElement
    if ismember(j,To_solve) || ismember(j,Active)
        if stopindic(j) == 0
            if iteracion == 1
                
                %% Phase-Field: Implicit method (LSCHEME)
                [Sol_Phasefieldtest,~,geometry{j}] = ...
                    LschemeSolution_phaseField(geometry{j},Sol_Phasefield_iter{j}.Pres,...
                    Sol_Phasefield_n{j}.Pres,react(j),in_tol_loop);
                
                if Sol_Phasefieldtest.radii < stopradii_loop
                    
                    Sol_Phasefield{j} = Stop_phasefield;
                    geometry1{j}      = Stop_geometry;
                    [Sol_Phasefield_n{j},~] = Projection_meshes1(Sol_Phasefield_n{j},geometry{j},geometry1{j});
                    
                    A_Efective(:,:,j) = Stop_A_Efective;
                    K_Efective(:,:,j) = Stop_K_Efective;
                    
                else
                    %% Refinement
                    [geometry1{j},geo_full,Sol_Phasefield_iter{j}] = AposterioriRefinement1...
                        (Sol_Phasefieldtest,Sol_Phasefield_iter{j},geometry{j});
                    
                    [geo_full]      = MFEM_preProcess_Micro(geo_full);
                    [geometry1{j}]  = MFEM_preProcess_Micro(geometry1{j});
                    
                    [Sol_Phasefield_n1,~] = Projection_meshes1(Sol_Phasefield_n{j},geometry{j},geo_full);
                    
                    %% Recalculate the phase field
                    
                    [Sol_Phasefield_full,~,~] = ...
                        LschemeSolution_phaseField(geo_full,...
                        Sol_Phasefield_iter{j}.Pres,Sol_Phasefield_n1.Pres,react(j),in_tol_loop);
                    
                    [Sol_Phasefield{j},~]      = Projection_meshes1(Sol_Phasefield_full,geo_full,geometry1{j});
                    [Sol_Phasefield_n{j},~]    = Projection_meshes1(Sol_Phasefield_n1,geo_full,geometry1{j});
                    %                     [Sol_Phasefield_iter{j},~] = Projection_meshes1(Sol_Phasefield_iter{j},geo_full,geometry1{j});
                    
                    %% CELL PROBLEMS -DIFFUSION & PERMEABILITY
                    [SolDiff1,SolDiff2,SolPerm1,SolPerm2] = ...
                        Solution_cellProblems(geometry1{j},Sol_Phasefield{j});
                    
                    %% Effective tensors -DIFFUSION & PERMEABILITY
                    [A_Efective(:,:,j), K_Efective(:,:,j)]= EfectiveTensor(geometry1{j},Sol_Phasefield{j}.Pres,...
                        SolDiff1,SolDiff2,SolPerm1,SolPerm2);
                end
                
            else
                % ITERACION >1 NO REFINEMENT
                [Sol_Phasefield{j},~,geometry1{j}] = ...
                    LschemeSolution_phaseField(geometry{j},Sol_Phasefield_iter{j}.Pres,...
                    Sol_Phasefield_n{j}.Pres,react(j),in_tol_loop);
                
                if Sol_Phasefield{j}.radii < stopradii_loop
                    Sol_Phasefield{j} = Stop_phasefield;
                    geometry1{j}      = Stop_geometry;
                    [Sol_Phasefield_n{j},~] = Projection_meshes1(Sol_Phasefield_n{j},geometry{j},geometry1{j});
                    
                    A_Efective(:,:,j) = Stop_A_Efective;
                    K_Efective(:,:,j) = Stop_K_Efective;
                else
                    %% CELL PROBLEMS -DIFFUSION & PERMEABILITY
                    [SolDiff1,SolDiff2,SolPerm1,SolPerm2] = ...
                        Solution_cellProblems(geometry1{j},Sol_Phasefield{j});
                    
                    %% Effective tensors -DIFFUSION & PERMEABILITY
                    [A_Efective(:,:,j), K_Efective(:,:,j)]= EfectiveTensor(geometry1{j},Sol_Phasefield{j}.Pres,...
                        SolDiff1,SolDiff2,SolPerm1,SolPerm2);
                end
            end
            % SI YA ANTES PARAMOS ENTONCES COPIAR
        else
            Sol_Phasefield{j} = Sol_Phasefield_iter{j};
            A_Efective(:,:,j) = A_Efective0(:,:,j);
            K_Efective(:,:,j) = K_Efective0(:,:,j);
            geometry1{j} = geometry{j};
        end
        Radii(j)    = Sol_Phasefield{j}.radii;
        Porosity(j) = Sol_Phasefield{j}.average;
    end
end

% Copy solution (DIRICHLET NODES)
if Dirichlet == 1
    for j=1:N_macro*2
        elem       = Micro_paste(j,1);
        elem_paste = Micro_paste(j,2:end);
        
        Radii(elem_paste)    = Radii(elem);
        Porosity(elem_paste) = Porosity(elem);
        stopindic(elem_paste,1)    = stopindic(elem,1);
        
        for ii = 1:(N_macro/2)-1
            Sol_Phasefield{elem_paste(ii),1} = Sol_Phasefield{elem};
            A_Efective(:,:,elem_paste(ii)) = A_Efective(:,:,elem);
            K_Efective(:,:,elem_paste(ii)) = K_Efective(:,:,elem);
            geometry1{elem_paste(ii),1}   = geometry1{elem};
        end
    end
else
    % COPY STRATEGY BETWEEN ACTIVE AND NON-ACTIVE
    if size(Active,1) < nElement
        for j=1:size(Active,1)
            if sum(simplycopy(:,j))>0
                elem       = Active(j);
                elem_paste = InActive(simplycopy(:,j)==1);
                
                Radii(elem_paste)    = Radii(elem);
                Porosity(elem_paste) = Porosity(elem);
                stopindic(elem_paste,1)    = stopindic(elem,1);
                
                for ii = 1:length(elem_paste)
                    Sol_Phasefield{elem_paste(ii),1} = Sol_Phasefield{elem};
                    A_Efective(:,:,elem_paste(ii))   = A_Efective(:,:,elem);
                    K_Efective(:,:,elem_paste(ii))   = K_Efective(:,:,elem);
                    geometry1{elem_paste(ii),1}      = geometry1{elem};
                    [Sol_Phasefield_n{elem_paste(ii)},~] = Projection_meshes1(Sol_Phasefield_n{elem_paste(ii)},...
                        geometry{elem_paste(ii)},geometry1{elem_paste(ii)});
                end
            end
        end
    end
end
end


