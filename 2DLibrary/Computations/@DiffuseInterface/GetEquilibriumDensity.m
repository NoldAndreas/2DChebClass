function [rho,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,rho,findTheta)
    Cn     = this.optsPhys.Cn;
    thetaEq = this.optsPhys.thetaEq;
    
    M      = this.IC.M;           
    Ind    = this.IC.Ind;    
    Diff   = this.IC.Diff;
    Int    = this.IC.Int;
    y2Max  = this.optsNum.PhysArea.y2Max;
    nParticles = this.optsPhys.nParticles;
    
    IntSubArea  = this.IntSubArea;          
%     Dy1 = this.optsNum.PhysArea.y2Max*cos(theta);
% 	bxArea          = struct('y1Min',-10+Dy1,'y1Max',10+Dy1,'N1',100,...
%                              'y2Min',this.optsNum.PhysArea.y2Max-5,...
%                              'y2Max',this.optsNum.PhysArea.y2Max,'N2',100);  
%     bxArea.N        = [100,100];
%     BX              = Box(bxArea);
%     IntBx           = BX.ComputeIntegrationVector();
%     IntSubArea = IntBx*this.IC.SubShapePts(BX.GetCartPts());
    
    
    if(isscalar(mu))
        mu = mu*ones(M,1);
    end                   
    
 %   y = fsolve(@f_eq,rho);
    
    %y        = NewtonMethod([0;rho],@f_eq);    %0;
 %   nParticles = -cos(theta)*(y2Max)^2;
    if((nargin > 4) && strcmp(findTheta,'findTheta'))
        y        = NewtonMethod([0;theta;rho],@f_eq_theta,0.001,100);    %0;

        muDelta  = y(1);
        theta    = y(2);
        rho      = y(3:end);
    else
        muDelta  = 0;
    end
    
  %  nParticles = -cos(theta)*(y2Max)^2;
    y        = NewtonMethod([muDelta;rho],@f_eq);    %0;
    muDelta  = y(1);    
    rho      = y(2:end);
    
    disp(['theta = ',num2str(theta*180/pi)]);
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
        rhoDiag = diag(rho_s);        
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - cos(thetaEq)*(1-rho_s(Ind.bottom)).^2/Cn;
        J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),...
                                Ind.normalBottom*Diff.grad-2*cos(thetaEq)*(1-rhoDiag(Ind.bottom,:))/Cn];
        
        %% Boundary condition for the density at y2max
        % 
        % $$(\cos(\theta),\sin \theta)^T \cdot {\nabla \rho} = 0$$
        %         
         E               = eye(M);
         ETop            = E(Ind.top,:);         
         topDirection    = [cos(theta)*ETop,sin(theta)*ETop];
         y(Ind.top)      = topDirection*(Diff.grad*rho_s); 
         J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad];                                  
%         topDirection    = [diag(this.IC.GetCartPts.y1_kv(Ind.top)),...
%                            diag(this.IC.GetCartPts.y2_kv(Ind.top))];
         %y(Ind.top)      = topDirection*(Diff.grad([Ind.top;Ind.top],:)*rho_s);         
         %J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad([Ind.top;Ind.top],:)];                
        
%         ys             = DoublewellPotential(rho_s,Cn) - Cn*Diff.Lap*rho_s;
%         y(Ind.top)     = - Diff.Lap(Ind.top,:)*((rho + rho_m).*mu_s) ...
%                          + Diff.div(Ind.top,:)*(repmat(ys,2,1).*(Diff.grad*rho_s)) + ... 
%                          + Cak*(zeta + 4/3*eta)*Diff.Lap(Ind.top,:)*Diff.div*uv;
%         J(Ind.top,:)   = - Diff.Lap(Ind.top,:)*diag(mu_s)
            
        %% Mass condition
        rhoM        = rho_s(Ind.left & Ind.top);
        rhoP        = rho_s(Ind.right & Ind.top);
        y           = [IntSubArea*(rho_s - (rhoM + rhoP)/2) - nParticles;y];
        
        Jint                      = IntSubArea;        
        Jint(Ind.left & Ind.top)  = IntSubArea(Ind.left & Ind.top)/2;
        Jint(Ind.right & Ind.top) = IntSubArea(Ind.right & Ind.top)/2;

        J           = [0,Jint;J];   
       %  J = J(:,2:end);
    end


    function [y,J] = f_eq_theta(x)
        
        dmu         = x(1);
        theta_s     = x(2);        
        mu_s        = mu + dmu;
        rho_s       = x(3:end); 
        [y,J]       = GetExcessChemPotential(rho_s,mu_s); 
        
        %% Boundary condition for the density at wall
        % 
        % $${\bf n}\cdot {\nabla \rho} = g$$
        %                                
        rhoDiag = diag(rho_s);        
        y(Ind.bottom)   = Ind.normalBottom*(Diff.grad*rho_s) - cos(thetaEq)*(1-rho_s(Ind.bottom)).^2/Cn;
        J(Ind.bottom,:) = [zeros(sum(Ind.bottom),1),...
                                Ind.normalBottom*Diff.grad-2*cos(thetaEq)*(1-rhoDiag(Ind.bottom,:))/Cn];
        
        %% Boundary condition for the density at y2max
        % 
        % $$(\cos(\theta),\sin \theta)^T \cdot {\nabla \rho} = 0$$
        %         
         E               = eye(M);
         ETop            = E(Ind.top,:);
         topDirection    = [cos(theta_s)*ETop,sin(theta_s)*ETop];
         topDirectionD   = [-sin(theta_s)*ETop,cos(theta_s)*ETop];
         y(Ind.top)      = topDirection*(Diff.grad*rho_s); 
         %J(Ind.top,:)    = [ topDirectionD*(Diff.grad*rho_s) ,topDirection*Diff.grad];                                  
         J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad]; 
         
         J_theta          = zeros(M,1);
         J_theta(Ind.top) = topDirectionD*(Diff.grad*rho_s);         
         J_inclTheta      = [J(:,1),J_theta,J(:,2:end)];
         
         
%         topDirection    = [diag(this.IC.GetCartPts.y1_kv(Ind.top)),...
%                            diag(this.IC.GetCartPts.y2_kv(Ind.top))];
         %y(Ind.top)      = topDirection*(Diff.grad([Ind.top;Ind.top],:)*rho_s);         
         %J(Ind.top,:)    = [zeros(sum(Ind.top),1),topDirection*Diff.grad([Ind.top;Ind.top],:)];                
        
%         ys             = DoublewellPotential(rho_s,Cn) - Cn*Diff.Lap*rho_s;
%         y(Ind.top)     = - Diff.Lap(Ind.top,:)*((rho + rho_m).*mu_s) ...
%                          + Diff.div(Ind.top,:)*(repmat(ys,2,1).*(Diff.grad*rho_s)) + ... 
%                          + Cak*(zeta + 4/3*eta)*Diff.Lap(Ind.top,:)*Diff.div*uv;
%         J(Ind.top,:)   = - Diff.Lap(Ind.top,:)*diag(mu_s)
            
        %% Mass condition       
        
        rhoM        = rho_s(Ind.left & Ind.top);
        rhoP        = rho_s(Ind.right & Ind.top);
        y           = [dmu;...
                       IntSubArea*(rho_s - (rhoM + rhoP)/2) - nParticles;...
                       y];
        
        Jint                      = IntSubArea;        
        Jint(Ind.left & Ind.top)  = IntSubArea(Ind.left & Ind.top)/2;
        Jint(Ind.right & Ind.top) = IntSubArea(Ind.right & Ind.top)/2;

        J           = [1,0,zeros(1,M);...
                       0,sin(theta_s)*(y2Max)^2,Jint;...
                       J_inclTheta];   
       %  J = J(:,2:end);
    end
  

end