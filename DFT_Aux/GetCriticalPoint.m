function [kBT_crit,rho_crit,mu_crit,p_crit] = GetCriticalPoint(optsPhys,initialGuess,output)

 %Load Data & Initializing    
    HS_f          = str2func(optsPhys.HSBulk);    	
    [h1,h2,alpha] = getV2(0,optsPhys.V2);    
    alpha
    
    if((nargin == 1) || isempty(initialGuess))
        initialGuess = [1;0.4;-3.4]; 
    end        
    
    %Solve Equation    
    opts = optimset('Display','off');
    [x,h1,flag] = fsolve(@f,initialGuess,opts);
    if(flag ~= 1)
        disp('Solution has not converged!!');
    end
    
    
    %Set Output
    kBT_crit = x(1);
	rho_crit = x(2); 
	mu_crit  = x(3); 
    
    [h1,CS_L] = HS_f(rho,kBT);
    p_crit    = -(kBT*rho_crit*(log(rho_crit) - 1) + alpha*rho_crit^2 + CS_L - mu_crit*rho_crit);
            
    if(((nargin < 3) || (output == true)))
        disp('***');
        disp('Values at Critical point:');
        disp(['kBT_crit = ',num2str(kBT_crit)]);
        disp(['rho_crit = ',num2str(rho_crit)]);
        disp(['mu_crit = ',num2str(mu_crit)]);    
        disp('***');
    end
    
    function z = f(x)
        kBT = x(1);
        rho = x(2); 
        mu  = x(3);         
        
        [mu_s,h1s,dmu_s,ddmu_s] = HS_f(rho,kBT);        
        
        dy   = kBT*(log(rho)) + 2*alpha*rho + mu_s - mu;
        ddy  = kBT/rho + 2*alpha + dmu_s;
        dddy = -kBT./(rho.^2) + ddmu_s;
        
        z = [dy;ddy;dddy];
    end    


end