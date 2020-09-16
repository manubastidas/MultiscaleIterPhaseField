%*********************************************************************
% PARAMETERS OF THE SIMULATION
%
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium


function getParameters()

%% Method
global  Dirichlet save_points N_macro Lstab

Dirichlet = 0;
N_macro   = 8;  %Partition Macro grid

Lstab     = 6;

save_points = [0.1, 0.1;0.5 0.25; 0.9 0.4];
 
%% GeneralParameters 
global  lambda gamma_par Ka n_param delta D_dif mu_f ustar 
% DO NOT TOUCH! 
lambda    = 0.08;
gamma_par = 1e-2;
D_dif   = 1;
mu_f    = 1;
Ka      = 25;
n_param = 10;
ustar   = 1;
delta   = 1E-4;

%% Time parameters 
global T_MAX Time 

T_MAX         = 0.25;
Time.tnSteps  = 25;
Time.dt       = T_MAX/Time.tnSteps;
Time.time_vec = 0+Time.dt:Time.dt:T_MAX+Time.dt;

Time.savetime = 1:Time.tnSteps; %unique([TRIN:TRIN:Time.tnSteps,Time.tnSteps]);


%% Inital conditions
global u_eq p_init
global u_down p_up p_down

u_eq   = 0.5;
p_init = 0;

% Boundary conditions
u_down = 0;
% u_up = 0;
p_up   = 0;
p_down = 0;

%% Initial micro mesh
global R N_micro fine_N Theta0 Theta1 Stop_R

% R = sqrt((1-solution.average)/pi)
R       = 0.4; % Initial ratio of the phasfields
Stop_R  = 0.1;
N_micro = 10;  %Partition micro grid

% NN -> CircularPhaseFieldT line10: compute intial pfield
fine_N  = 200; 

Theta0  = lambda/3; %Necesary resolution
Theta1  = 2; %Micro-scale refinement parameter

%% Mesh refinement
global in_tol out_tol macro_tolr macro_lambda

in_tol  = 1E-8; %Tolerance non-linear solver
out_tol = 1E-6; %Tolerance multi-scale (use 1 for the micro experiments)

macro_lambda = 0.1;
macro_tolr   = 0; %  0 = no macro refinement

save('Parameters.mat')