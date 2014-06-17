function rho = GetEquilibriumDensity(this,mu,theta,nParticles,rho)
    Cn     = this.optsPhys.Cn;
    g      = this.optsPhys.g;    
    
    M      = this.IC.M;           
    Ind    = this.IC.Ind;    
    Diff   = this.IC.Diff;
    Int    = this.IC.Int;
    
    if(isscalar(mu))
        mu = mu*ones(M,1);
    end                   
    
    y     = NewtonMethod([0;rho],@f_eq);
    rho   = y(2:end);
    
	
    function [mu_s,J] = GetExcessChemPotential(rho_s,mu_offset)    
        [dW,~,ddW]    = DoublewellPotential(rho_s,Cn);
        mu_s          = dW - Cn*(Diff.Lap*rho_s) - mu_offset;                                   
        J             = diag(ddW) - Cn*Diff.Lap;   
        J             = [-ones(length(rho_s),1),J];
    end
    function [y,J] = f_eq(x)
        
        dmu         = x(1);
        rho_s       = x(2:end);
        [y,J]       = GetExcessChemPotential(rho_s,mu+dmu); 
        
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
        y           = [Int*rho_s - nParticles;y];
        J           = [0,Int;J];   
    end
  

end