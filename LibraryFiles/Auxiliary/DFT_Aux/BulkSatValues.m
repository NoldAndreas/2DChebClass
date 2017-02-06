function [rhoGas_sat,rhoLiq_sat,mu_sat,p] = BulkSatValues(optsPhys,intitialGuess,output)
% optsPhys - struct with 'V2.V2DV2','HSBulk','kBT'

    %Load Data & Initializing
    if((~isfield(optsPhys,'V2')) || ...
            (isfield(optsPhys.V2,'V2DV2') && strcmp(optsPhys.V2.V2DV2,'zeroPotential')))
        alpha = 0;
    else   
        [h1,h2,alpha] = getV2(0,optsPhys.V2);            
    end
        
    HS_f  = str2func(optsPhys.HSBulk);
    kBT   = optsPhys.kBT;
    
    if(nargin == 1)
        intitialGuess = [0.02;0.85;-2.0];
        %intitialGuess = [-0.02i;0.85;-2.0];
    end        
    
    %Solve Equation    
    opts = optimset('Display','off');
    x = fsolve(@fBulk,intitialGuess,opts);
    
    
    %Set Output
    rhoGas_sat = x(1);
    rhoLiq_sat = x(2);
    mu_sat     = x(3);
    
    [h1,CS_L] = HS_f(rhoLiq_sat,kBT);
    p         = -(kBT*rhoLiq_sat*(log(rhoLiq_sat) - 1) + alpha*rhoLiq_sat^2 + CS_L - mu_sat*rhoLiq_sat);
    
    if(abs(rhoGas_sat-rhoLiq_sat)<1e-8)
        cprintf('*r','BulkSatValues: No Liquid and Gas Phas found. Only one Phase!');
    end
            
    if((nargin < 3) || output)
        disp('***');
        disp('Saturation Values:');
        disp(['mu = ',num2str(mu_sat)]);
        disp(['rho_gas = ',num2str(rhoGas_sat)]);
        disp(['rho_liq = ',num2str(rhoLiq_sat)]);
        disp(['p_Bulk  = ',num2str(p)]);
        disp('***');
    end
    
    function z = fBulk(x)
        rho1        = x(1); 
        rho2        = x(2); 
        mu          = x(3);
        
        [muCS1,CS1] = HS_f(rho1,kBT);
        [muCS2,CS2] = HS_f(rho2,kBT);
        
        dy1 = kBT*(log(rho1)) + 2*alpha*rho1 + muCS1 - mu;
        dy2 = kBT*(log(rho2)) + 2*alpha*rho2 + muCS2 - mu;
        y1  = kBT*rho1*(log(rho1) - 1) + alpha*rho1^2 + CS1 - mu*rho1;
        y2  = kBT*rho2*(log(rho2) - 1) + alpha*rho2^2 + CS2 - mu*rho2;
        
        y1My2 = y1 - y2;
        
        z = [dy1;dy2;y1My2];
    end    

end