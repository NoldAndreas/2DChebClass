function rho = GetEquilibriumDensity(this,mu,theta,nParticles,rho)
    Cn     = this.optsPhys.Cn;
    g      = this.optsPhys.g;
    rhoInf = -1;
    
    M      = this.IC.M;           
    Ind    = this.IC.Ind;    
    Diff   = this.IC.Diff;
    Int    = this.IC.Int;
         
    bulkSolve = (~Ind.right & ~Ind.left);
    
    if(isscalar(mu))
        mu = mu*ones(M,1);
    end
            
    rhoInf = NewtonMethod(rhoInf,@f_eq_inf);

    y     = NewtonMethod([0;rho(bulkSolve)],@f_eq);                 
    rho   = GetFullRho(y(2:end));
    
	function [y,dy] = f_eq_inf(rhoInf)
        muInf       = mu(Ind.right);                
        [~,y,dy]    = DoublewellPotential(rhoInf,Cn);
        y           = DoublewellPotential(rhoInf,Cn) - muInf(1);        
    end
    function [mu_s,J] = GetExcessChemPotential(rho_s,mu_offset)    
        [dW,W,ddW]    = DoublewellPotential(rho_s,Cn);
        mu_s          = dW - Cn*(Diff.Lap*rho_s) - mu_offset;                                   
        J             = diag(ddW) - Cn*Diff.Lap;
        J             = [-ones(length(rho),1),J];  
    end
    function [y,J] = f_eq(x)
        %solves for T*log*rho + Vext          
        mu_s        = x(1);        
        rho_s       = GetFullRho(x(2:end));
                
        [y,J] = GetExcessChemPotential(rho_s,mu+mu_s);        
        
        
        %% Boundary condition for the density at wall
        % 
        % $${\bf n}\cdot {\nabla \rho} = g$$
        %         
        
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - g;
        J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),Ind.normalBottom*Diff.grad];
        
        %% Boundary condition for the density at y2max
        % 
        % $$(\cos(\theta),\sin \theta)^T \cdot {\nabla \rho} = 0$$
        % 
        E               = eye(M);
        ETop            = E(Ind.top,:);
        topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
        y(Ind.top)      = topDirection*(Diff.grad*rho_s); %Ind.normalTop*(Diff.grad*rho_s); %TODO!!!
        J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad];
        
        %% Mass condition
        y           = [Int*rho_s - nParticles;y(bulkSolve)];
        J           = [0,Int(bulkSolve);J(bulkSolve,[true;bulkSolve])];
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        %y           = y(bulkSolve);
        %J           = J(bulkSolve,[false;bulkSolve]);
        %y           = [Int_SubOnFull*rho_s - nParticles;y(bulkSolve)];
        %J           = [Int_SubOnFull(bulkSolve);J(bulkSolve,[false;bulkSolve])];
        
    end

    function rho_Full = GetFullRho(rho_s)
        rho_Full             = zeros(M,1);
        rho_Full(bulkSolve)  = rho_s;
        rho_Full(Ind.right)  = rhoInf;
        rho_Full(Ind.left)   = 1; 
    end    

end