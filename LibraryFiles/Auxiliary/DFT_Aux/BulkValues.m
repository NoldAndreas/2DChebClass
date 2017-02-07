function [rhoGas_eq,rhoLiq_eq,pLiq,pGas] = BulkValues(mu,optsPhys,intitialGuess,output)
    
    

    %Load Data & Initializing
    if((~isfield(optsPhys,'V2')) || strcmp(optsPhys.V2.V2DV2,'zeroPotential'))
        alpha = 0;
    else
        [h_1,h_2,alpha] = getV2(0,optsPhys.V2);        
    end
    
    %Phi_r = str2func(optsPhys.V2DV2);        
    %[h1,h2,alpha] = Phi_r(0);    
    
    HS_f  = str2func(optsPhys.HSBulk);
    kBT   = optsPhys.kBT;    
    
    if((nargin < 3) || (isempty(intitialGuess)))
        intitialGuess = [0.01;1.1];
    end
    
    %Solve Equation    
    opts = optimset('Display','off','TolFun',1e-14,'TolX',1e-14);

    rhoGas_eq = fsolve(@fBulk,intitialGuess(1),opts);
    rhoLiq_eq = fsolve(@fBulk,intitialGuess(2),opts);   
    
    [h1,CS_L] = HS_f(rhoLiq_eq,kBT);
    [h1,CS_G] = HS_f(rhoGas_eq,kBT);
    pLiq = -(kBT*rhoLiq_eq*(log(rhoLiq_eq) - 1) + alpha*rhoLiq_eq^2 + CS_L - mu*rhoLiq_eq);
    pGas = -(kBT*rhoGas_eq*(log(rhoGas_eq) - 1) + alpha*rhoGas_eq^2 + CS_G - mu*rhoGas_eq);

    if((nargin < 4) || output)
        disp('***');
        disp('Equilibrium Values (for given mu):');
        disp(['mu = ',num2str(mu)]);
        disp(['rho_gas = ',num2str(rhoGas_eq)]);
        disp(['rho_liq = ',num2str(rhoLiq_eq)]);    
        disp(['p_liq = ',num2str(pLiq)]);    
        disp(['p_gas = ',num2str(pGas)]);
        disp('***');
    end
%    disp(['Equilibrium Values at mu = ',num2str(mu),' (rhoLiq,rhoGas) = (',...
%        num2str(rhoLiq_eq),' , ' , num2str(rhoGas_eq),')']);   
    
    
    function dy = fBulk(rho)
        dy = kBT*(log(rho)) + 2*alpha*rho + HS_f(rho,kBT) - mu;        
    end
    
end