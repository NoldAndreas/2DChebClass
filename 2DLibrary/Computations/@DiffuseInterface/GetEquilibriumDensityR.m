function [rho,muDelta] = GetEquilibriumDensityR(this,mu,theta,nParticles,rho,ptC)
    Cn     = this.optsPhys.Cn;
    g      = this.optsPhys.g;    
    
    M      = this.IC.M;           
    Ind    = this.IC.Ind;    
    Diff   = this.IC.Diff;
    Int    = this.IC.Int;
    PtsCart = this.IC.GetCartPts;
    IntSubArea  = this.IntSubArea;   
    
    if(isscalar(mu))
        mu = mu*ones(M,1);
    end                   
    
 %   y = fsolve(@f_eq,rho);
    
    y        = NewtonMethod([0;rho],@f_eq);    %0;
    
    muDelta  = y(1);
    rho      = y(2:end);
    
    disp(['Delta mu = ',num2str(muDelta)]);
    
    %
	
    function [mu_s,J] = GetExcessChemPotential(rho_s,mu_offset)    
        [dW,~,ddW]    = DoublewellPotential(rho_s,Cn);
        mu_s          = dW - Cn*(Diff.Lap*rho_s) - mu_offset;                                   
        J             = diag(ddW) - Cn*Diff.Lap;   
        J             = [-ones(length(rho_s),1),J];
    end
    function [y,J] = f_eq(x)
        
        dmu         = x(1);
        mu_s        = mu + dmu;
        rho_s       = x(2:end);         
        [y,J]       = GetExcessChemPotential(rho_s,mu_s); 
        
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
         ETop            = eye(sum(Ind.top));         
         topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
         
         y1CT  = PtsCart.y1_kv(Ind.top);                  
         y2CT  = PtsCart.y2_kv(Ind.top);
         
         mergeF          = y1CT.^2./(10000 + y1CT.^2);
         mergeF(1) = 1; mergeF(end) = 1;
         mergeF          = [diag(mergeF),diag(mergeF)];
        
         y1CT = y1CT - ptC.y1;
         y2CT = y2CT - ptC.y2;
         
         length          = sqrt(y1CT.^2+y2CT.^2);
         length(1) = 1; length(end) = 1;
         y1CT(1) = 1; y1CT(end) = 1;
         curvDirection   = [-diag(y2CT./length),diag(y1CT./length)];
                  
         direction       = mergeF.*topDirection + (1-mergeF).*curvDirection;         
         
         y(Ind.top)      = direction*(Diff.grad([Ind.top;Ind.top],:)*rho_s); 
         J(Ind.top,:)    = [zeros(sum(Ind.top),1),direction*Diff.grad([Ind.top;Ind.top],:)];
            
        %% Mass condition
        rhoM        = rho_s(Ind.left & Ind.top);
        rhoP        = rho_s(Ind.right & Ind.top);
        y           = [IntSubArea*(rho_s - (rhoM + rhoP)/2) - nParticles;y];
        
        Jint                      = IntSubArea;        
        Jint(Ind.left & Ind.top)  = IntSubArea(Ind.left & Ind.top)/2;
        Jint(Ind.right & Ind.top) = IntSubArea(Ind.right & Ind.top)/2;

        J           = [0,Jint;J];   
         %J = J(:,2:end);
    end
  

end