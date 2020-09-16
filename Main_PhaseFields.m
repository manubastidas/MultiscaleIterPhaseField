%*********************************************************************
% Main Code for Adaptime multi-scale iterative solver 
% for a two-scale phase-field model 
%
%*********************************************************************
%
%***------------------------------------
%***Inputs: Parameters phase field
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

clear
close all
clc

tic;
getParameters()

global N_macro Dirichlet macro_tolr

% Write name of the results
% name = sprintf('Example1_Lcoup1_L6');

CurrentDir = pwd();
% Write name of the folder
src_name = 'Run_all';
if ~strcmp(CurrentDir(end-length(src_name):end),['\',src_name]) && ~strcmp(CurrentDir(end-length(src_name):end), ['/',src_name])
    error('You have to be in the folder src to run the code!');
else
    root_folder = CurrentDir(1:end-length(src_name));
end
resultDir = [root_folder,name];

mkdir(resultDir)

parfile = fullfile(resultDir,['Param_',name,'.mat']);
save(parfile)

global Time out_tol Stop_R u_down save_points R

%% Macro Geometry
% Basic geometry [0,1]\times[0,1/2]
if N_macro>1
    %  SOLVING MORE THAN ONE MACRO POINT
    [xx,yy]             = meshgrid(0:1/N_macro:1,0:1/N_macro:1/2);
    Macrogeo            = struct();
    Macrogeo.coordinate = [xx(:),yy(:)];
    Macrogeo.element    = delaunay(xx,yy);
    Macrogeo.TR         = triangulation(Macrogeo.element,Macrogeo.coordinate);
    
    % Macrogeo features
    [Macrogeo]           = myedge(Macrogeo);
    [MatC,MatD,Macrogeo] = MFEM_preProcess_Macro(Macrogeo,N_macro);
else
    Macrogeo.nElement = 1;
    Macrogeo.nnodes   = 1;
    Macrogeo.area     = 1;
end
% Dirichlet conditions: Only solve one layer of micro scale
if Dirichlet ==1
    % For example projection 1D
    [Macrogeo,Micro_paste] = DirichletConditions(Macrogeo);
else
    Macrogeo.uDir_Gamma = [];
    Macrogeo.pDir_Gamma = [];
    Micro_paste         = (1:Macrogeo.nElement)';
end

field = sprintf('time%i',0);
[MacroSol.(field)] = InitialSolution(Macrogeo);

%% Micro scale
A_Efective = struct();
K_Efective = struct();
A_Efective.(field) = ones(2,2,Macrogeo.nElement);
K_Efective.(field) = ones(2,2,Macrogeo.nElement);

Sol_Phasefield = struct();

% Basic Initial phase-field (Can be change for squares or other geometries)
[Microgeometry1,Sol_Phasefield.(field),Radii.(field),Porosity.(field),...
    Stop_geometry,Stop_Phasefield,Stop_Radii,Stop_Porosity,Stop_A_Efective,Stop_K_Efective] =...
    InitialPhaseField(Macrogeo,R,Stop_R);
% disp(Stop_Radii)
% Indicator of the rattio less than R0
stopindic    = zeros(Macrogeo.nElement,1);

% FOR THE Macro Refinement
% VectSimilarities = [MacroSol.(field).u,Porosity.(field)'];
Active   = cell(Time.tnSteps,1);
InActive = cell(Time.tnSteps,1);

Micro_dof = zeros(Time.tnSteps,1);

simplycopy  = cell(Time.tnSteps,1);
Distance = zeros(Macrogeo.nElement);
if Dirichlet ==1
    Active{1}   = Micro_paste(:,1) ;
elseif macro_tolr >0
    Active{1}   = randi(1:Macrogeo.nElement);
    [Active{1},InActive{1},Distance,simplycopy{1},Micro_dof(1,:),stopindic] = ...
        Macro_adaptivity(Sol_Phasefield.(field),Microgeometry1,...
        MacroSol.(field).u,Distance,Active{1},Macrogeo,Porosity.(field),stopindic);
else
    Active{1} = (1:Macrogeo.nElement)';
    InActive{1}  = [];
    simplycopy{1} = zeros();
end 
    
% Auxiliar variables
Microgeometry         = struct();
Microgeometry.(field) = Microgeometry1;
Sol_Phasefield0       = Sol_Phasefield.(field);
% Microgeometry0        = Microgeometry1;

MicrogeometrySave     = struct();
Sol_PhasefieldSave    = struct();

field0  = field;

% Errors
MacroError_iter = cell(Time.tnSteps,1);
PhaseError      = zeros(Time.tnSteps,Macrogeo.nElement,0);
PhaseError_iter = cell(Time.tnSteps,1);
errorIter       = cell(Time.tnSteps,1);

init_t = toc;

% TO REPORT
it        = ones(Time.tnSteps,1);
time_t = zeros(Time.tnSteps,1);

%
fprintf('\n Time %i/%i \n',0,Time.tnSteps)
for kk=1 %:Time.tnSteps
    tic
    
    field  = sprintf('time%i',kk);
    
    %     it = 1;
    errorIter{kk}(1) = inf;
    
    %     Microgeometry1        = Microgeometry.(field0);
    Sol_Phasefield_iter   = Sol_Phasefield.(field0);
    Sol_Phasefield_n      = Sol_Phasefield.(field0);
    A_Efective_iter       = A_Efective.(field0);
    K_Efective_iter       = K_Efective.(field0);
    MacroSol_iter         = MacroSol.(field0);
    Porosity_iter         = Porosity.(field0);
    
    PhaseError_iter{kk}(it(kk)) = inf;
    
    while errorIter{kk}(it(kk)) > out_tol
        
        Porosity_aux = Porosity_iter;
        
        [Radii_iter,Porosity_iter,Sol_Phasefield_iter,Sol_Phasefield_n,...
            A_Efective_iter,K_Efective_iter,Microgeometry1] = ...
            MicroCellProblems_Lscheme(Microgeometry1,Macrogeo,MacroSol_iter.react,...
            Sol_Phasefield_iter,Sol_Phasefield_n,...
            A_Efective_iter,K_Efective_iter,stopindic,...
            Stop_Phasefield,Stop_geometry,Stop_A_Efective,...
            Stop_K_Efective,it(kk),Micro_paste,Active{kk},InActive{kk},simplycopy{kk});
        
        aux_u = MacroSol_iter;
        
        if N_macro>1
            % MACRO PROBLEM 2 -- (for p)
            [MacroSol_iter] = MacroProb2(Macrogeo,K_Efective_iter,MatC);
            
            % MACRO PROBLEM 1 -- (for u)
            [MacroSol_iter] = MacroProb1(MacroSol_iter,Macrogeo,A_Efective_iter,...
                Porosity.(field0),Porosity_iter,...
                MacroSol.(field0).u,aux_u.u,MatC,MatD);
        else
            MacroSol_iter.u          = u_down;
            MacroSol_iter.react(:,1) = (MacroSol_iter.u.^2)/(0.5^2)-1;
        end
        
        it(kk) = it(kk)+1;
        %         PhaseError(kk,:,it)     = sqrt(phaseerror)';
        
        MacroError_iter{kk}(it(kk)) = sqrt(sum((aux_u.u-MacroSol_iter.u).^2.*Macrogeo.area'));
        PhaseError_iter{kk}(it(kk)) = sqrt(sum((Porosity_aux-Porosity_iter).^2.*Macrogeo.area));
        errorIter{kk}(it(kk))       = PhaseError_iter{kk}(it(kk));
        
        % % --------------------------------------------
        % Report - error
                fprintf('\n Time %i - Iter %i - %1.2E',kk,it(kk),errorIter{kk}(it(kk)))
        %         disp([PhaseError_iter{kk}(it(kk)),MacroError_iter{kk}(it(kk)),PhaseError_iter{kk}(it(kk))+MacroError_iter{kk}(it(kk))])
        
    end
    
    A_Efective.(field)     = A_Efective_iter;
    K_Efective.(field)     = K_Efective_iter;
    Porosity.(field)       = Porosity_iter;
    Radii.(field)          = Radii_iter;
    MacroSol.(field)       = MacroSol_iter;
    
    Microgeometry.(field)  = Microgeometry1;
    Sol_Phasefield.(field) = Sol_Phasefield_iter;
    
    if ismember(kk,Time.savetime)
        if N_macro>1
            for ss = 1:size(save_points,1)
                [~,one] = min(sum((Macrogeo.bari-repmat(save_points(ss,:),Macrogeo.nElement,1)).^2,2));
                MicrogeometrySave.(field){ss}  = Microgeometry.(field){one};
                Sol_PhasefieldSave.(field){ss} = Sol_Phasefield.(field){one};
            end
        else
            MicrogeometrySave.(field){1}  = Microgeometry.(field){1};
            Sol_PhasefieldSave.(field){1} = Sol_Phasefield.(field){1};
        end
    end
    
    %% Macro refinement
    if Dirichlet == 1 || N_macro==1 || macro_tolr ==0
        Active{kk+1} = Active{kk};
        if N_macro == 1
            Micro_dof(kk+1)= Microgeometry.(field){1}.nElement;
        elseif Dirichlet == 1 || macro_tolr ==0
            Micro_dof(kk+1)= 0;
            for rs = 1:length(Active{kk})
                Micro_dof(kk+1)= Micro_dof(kk+1) + Microgeometry.(field){Active{kk}(rs)}.nElement;
            end
        end
    else
        %         VectSimilarities = [MacroSol.(field).u,Porosity.(field)'];
        [Active{kk+1},InActive{kk+1},Distance,simplycopy{kk+1},Micro_dof(kk+1,:),stopindic] = ...
            Macro_adaptivity(Sol_Phasefield.(field),Microgeometry.(field),...
            MacroSol.(field).u,Distance,Active{kk},Macrogeo,Porosity.(field),stopindic);
    end
    
    field0 = field;
    fprintf('\n Time %i/%i \n',kk,Time.tnSteps)
    
    % Time to report
    time_t(kk) = toc;
end

load('Parameters.mat')
matfile = fullfile(resultDir, ['Results_',name,'.mat']);
save(matfile)

clc
disp(' --------------- !!!')
disp(name)
disp(' --------------- !!!')

% run report_make