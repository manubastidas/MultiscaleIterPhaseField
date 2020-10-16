%% source code

function [fplus,fminus,positions] = source_Lscheme(phi,phi_n,fu)

load('Parameters.mat','lambda','gamma_par')

x1 = -((-8*fu*lambda - 96*gamma_par ) + sqrt((8*fu*lambda + 96*gamma_par )^2 + 384*gamma_par*(-4*fu*lambda - 16*gamma_par )))/(192*gamma_par );
x2 = -((-8*fu*lambda - 96*gamma_par ) - sqrt((8*fu*lambda + 96*gamma_par )^2 + 384*gamma_par*(-4*fu*lambda - 16*gamma_par )))/(192*gamma_par );

positions = phi>=x1&phi<=x2; %% POSITIVE DERIVATIVE

split = 2;

if split == 1
    
    fprev   = -4*phi.*(1-phi).*(4*gamma_par*(1-2*phi)+ fu*lambda);
    fprev_n = -4*phi_n.*(1-phi_n).*(4*gamma_par*(1-2*phi_n)+ fu*lambda);
    
    fplus            = zeros(size(fprev_n));
    fplus(positions) = fprev(positions);
    fminus           = fprev;
    fminus(positions)= 0;
else
    
    %     funcion   = @(phi) 16*gamma*(1-phi).*phi.*(2*phi-1)-lambda*4*phi.*(1-phi)*fu(i);
    derivada  = @(phi)  16*gamma_par*(-1 + 6*phi - 6*phi.^2) + 2*lambda*(4*fu.*phi-2*fu);
    %     derivada1 = @(phi)  -96*gamma*(2*phi-1) - 8*lambda;
    
    derivadaplus  = @(phi) max(derivada(phi),0);
    derivadaminus = @(phi) min(derivada(phi),0);
    
    Fplus  = @(phi)integral(derivadaplus,0,phi);
    Fminus = @(phi)integral(derivadaminus,0,phi);
    
    fplus  = arrayfun(Fplus,phi_n);
    fminus = arrayfun(Fminus,phi);
end

