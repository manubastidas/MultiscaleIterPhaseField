%*********************************************************************
% Splitting and evaluation of the non-linear term F
%*********************************************************************
%
%***------------------------------------
% Manuela Bastidas - 2020
% Hasselt University, Belgium

function [ff,x1,x2] = source_Lscheme(phi,fu)

load('Parameters.mat','lambda','gamma_par')
% global lambda gamma  

% gamma = 1e-2;
% ustar = 1;
% f = @(x) (x./0.5).^2-1;

% fu = f(u);
ff = -4*phi.*(1-phi).*(4*gamma_par*(1-2*phi)+ fu*lambda);
x1 = -((-8*fu*lambda - 96*gamma_par ) + sqrt((8*fu*lambda + 96*gamma_par )^2 + 384*gamma_par*(-4*fu*lambda - 16*gamma_par )))/(192*gamma_par );
x2 = -((-8*fu*lambda - 96*gamma_par ) - sqrt((8*fu*lambda + 96*gamma_par )^2 + 384*gamma_par*(-4*fu*lambda - 16*gamma_par )))/(192*gamma_par );

% fp = 4*f(u)*lambda*(2*phi-1) - 16*gamma*(1 - 6*phi + 6*phi.^2);

%% Analisis of fp
% 
% fu = f(0);
% funcion = @(phi) -4*phi.*(1-phi).*(4*gamma*(1-2*phi)+ fu*lambda); 
% derivada = @(phi) 4*fu*lambda*(2*phi-1) - 16*gamma*(1 - 6*phi + 6*phi.^2);
% % derivada2 = @(phi) 8*fu*lambda - 16*gamma*(12*phi -6);
% derivada2 = @(phi) (8*fu*lambda + 96*gamma) - 192*gamma*phi;
% fplot(funcion,[0,1])
% hold on
% fplot(derivada,[0,1],'r')
% hold on
% fplot(derivada2,[0,1],'g')
% 
% % puntos donde cambia la derivada
% x = -((-8*fu*lambda - 96*gamma ) + sqrt((8*fu*lambda + 96*gamma )^2 + 384*gamma*(-4*fu*lambda - 16*gamma )))/(192*gamma )
% plot(x,0,'o')
% x2 = -((-8*fu*lambda - 96*gamma ) - sqrt((8*fu*lambda + 96*gamma )^2 + 384*gamma*(-4*fu*lambda - 16*gamma )))/(192*gamma )
% plot(x2,0,'o')
% grid on
% 
% %% choose L
% u = linspace(0,1);
% fu = f(u)';
% 
% posmaximo = [zeros(100,1) (96*gamma + 8*fu*lambda)./(192*gamma) ones(100,1)];
% der = [derivada(posmaximo(:,1)) derivada(posmaximo(:,2)) derivada(posmaximo(:,3))]
% maximo = max(der,2);
% 
% plot(maximo,derivada2(maximo),'o')