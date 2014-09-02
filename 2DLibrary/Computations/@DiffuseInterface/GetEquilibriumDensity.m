function [rho,theta,muDelta] = GetEquilibriumDensity(this,mu,theta,rho,findTheta)
    Cn     = this.optsPhys.Cn;
    thetaEq = this.optsPhys.thetaEq;
    if(length(thetaEq) == 1)
        thetaEq = thetaEq*[1,1];
    end
    
    M      = this.IC.M;           
    Ind    = this.IC.Ind;    
    Diff   = this.IC.Diff;        
    nParticles = this.optsPhys.nParticles;            
    
    IntSubArea  = this.IntSubArea;          
            
    if(isscalar(mu))
        mu = mu*ones(M,1);
    end                   
        
    a_direction = Set_a_direction(theta);
    
    if((nargin > 4) && strcmp(findTheta,'findTheta'))
        y        = NewtonMethod([theta;0;rho],@f_eq_theta,1e-6,100,0.5);    %0;

        theta    = y(1);
        muDelta  = y(2);        
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
	
    
    function [a_direction,a_direction_theta] = Set_a_direction(theta)
        a_y1             = zeros(M,1); a_y2             = zeros(M,1);               
        
        a_y1(Ind.bottom) = 0;          
        a_y2(Ind.bottom) = -1;           %BC at wall:        
        
        a_y1(Ind.top) = 0;          
        a_y2(Ind.top) = 1;
        
        a_y1(Ind.fluidInterface)    = cos(theta); 
        a_y2(Ind.fluidInterface)    = sin(theta);  %BC at fluid interface
        a_direction      = [diag(a_y1),diag(a_y2)];    

        a_y1  = zeros(M,1);
        a_y2  = zeros(M,1);
        a_y1(Ind.fluidInterface)    = -sin(theta);
        a_y2(Ind.fluidInterface)    = cos(theta);  %BC at fluid interface        
        a_direction_theta = [diag(a_y1),diag(a_y2)];
        
    end    
    function [g,Jg] = get_g(rho_s)
        g   = zeros(M,1);
        Jg  = zeros(M,M);
        rhoDiag         = diag(rho_s);
        
        g(Ind.bottom)    = cos(thetaEq(1))*(1-rho_s(Ind.bottom)).^2/Cn;                            
        Jg(Ind.bottom,:) = 2*cos(thetaEq(1))*(1-rhoDiag(Ind.bottom,:))/Cn;
        
        g(Ind.top)       = cos(thetaEq(2))*(1-rho_s(Ind.top)).^2/Cn;                            
        Jg(Ind.top,:)    = 2*cos(thetaEq(2))*(1-rhoDiag(Ind.top,:))/Cn;
        
        g(Ind.fluidInterface,:)  = 0;
        Jg(Ind.fluidInterface,:) = 0;
    end

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
        
        [g,Jg] = get_g(rho_s);
        
        % Boundary condition for the density
        topBottom       = Ind.bottom | Ind.top;
        y(topBottom)   = a_direction(topBottom,:)*(Diff.grad*rho_s) - g(topBottom);
        J(topBottom,:) = [zeros(sum(topBottom),1),...
                                a_direction(topBottom,:)*Diff.grad-Jg(topBottom,:)];
                           
        % Mass condition
        rhoM        = rho_s(Ind.left & Ind.top);
        rhoP        = rho_s(Ind.right & Ind.top);
        y           = [IntSubArea*(rho_s - (rhoM + rhoP)/2) - nParticles;y];
        
        Jint                      = IntSubArea;        
        Jint(Ind.left & Ind.top)  = IntSubArea(Ind.left & Ind.top)/2;
        Jint(Ind.right & Ind.top) = IntSubArea(Ind.right & Ind.top)/2;

        J           = [0,Jint;J];   
    end

    function [y,J] = f_eq_theta(x)
        
        theta_s     = x(1);  
        dmu         = x(2);                    
        rho_s       = x(3:end); 
        
        [y1,J1]     = f_eq(x(2:end));
        y           = [dmu;y1];
                
        [a_direction,a_direction_theta] = Set_a_direction(theta_s);        
        
        J_theta          = zeros(M,1);
        J_theta(Ind.top) = a_direction_theta(Ind.top,:)*(Diff.grad*rho_s);         
        
        J                = [0,1,zeros(1,M);...
                            0,J1(1,:);...
                            J_theta,J1(2:end,:)];
         
    end
  

end