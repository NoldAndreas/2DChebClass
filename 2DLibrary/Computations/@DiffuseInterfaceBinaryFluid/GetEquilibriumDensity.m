function [phi,theta] = GetEquilibriumDensity(this,theta,phi)%,findTheta)

    Cn      = this.optsPhys.Cn;
    m       = this.optsPhys.mobility;
    thetaEq = this.optsPhys.thetaEq;
    if(length(thetaEq) == 1)
        thetaEq = thetaEq*[1,1];
    end
    if(nargin < 2)
        theta = pi/2;
    end
    if(nargin < 3)
        phi = InitialGuessRho(this);
    end
    
    M       = this.IC.M;           
    Ind     = this.IC.Ind;    
    Diff    = this.IC.Diff;                    
    uv      = zeros(2*M,1);
        
    a_direction = Set_a_direction(theta);
    
    Convect    = [diag(uv(1:end/2)),diag(uv(1+end/2:end))]*Diff.grad;
    
    phi        = NewtonMethod(phi,@f_eq);    %0;
    %phi        = fsolve(@f_eq,phi);    %0;
    disp(['theta = ',num2str(theta*180/pi)]);    
	this.rho = phi;
    
    function [a_direction,a_direction_theta] = Set_a_direction(theta)
        a_y1             = zeros(M,1); 
        a_y2             = zeros(M,1);               
        
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
        
        g(Ind.bottom)    = cos(thetaEq(1))*(1-(rho_s(Ind.bottom)).^2)/Cn;                            
        Jg(Ind.bottom,:) = - 2*cos(thetaEq(1))*rhoDiag(Ind.bottom,:)/Cn;
        
        g(Ind.top)       = cos(thetaEq(2))*(rho_s(Ind.top)).^2/Cn;                            
        Jg(Ind.top,:)    = - 2*cos(thetaEq(2))*(1-rhoDiag(Ind.top,:))/Cn;
        
        g(Ind.fluidInterface)    = 0;
        Jg(Ind.fluidInterface,:) = 0;
    end

    function [mu_s,J] = GetExcessChemPotential(rho_s)
        [dW,~,ddW]    = DoublewellPotential(rho_s,Cn);
        mu_s          = dW - Cn*(Diff.Lap*rho_s);  
        J             = diag(ddW) - Cn*Diff.Lap;   
        J             = [-ones(length(rho_s),1),J];
    end

    function [y,J] = f_eq(phi_s)                       
        
        %y      = m*Diff.Lap*GetExcessChemPotential(phi_s) - Convect*phi_s;
        absGradRho  = ( (Diff.Dy1*phi_s).^2 + (Diff.Dy2*phi_s).^2 );
        
        y      = - Convect*phi_s ...
                 + 12*m/Cn*phi_s.*absGradRho ...
                 + 2*m/Cn*(3*phi_s.^2 -1).*(Diff.Lap*phi_s)...
                 - m*Cn*Diff.Lap2*phi_s;
        
        J      = 12*m/Cn*diag(absGradRho + phi_s.*(Diff.Lap*phi_s))...
                  - Convect ...
                  + 24*m/Cn*diag(phi_s)*(diag(Diff.Dy1*phi_s)*Diff.Dy1 + diag(Diff.Dy2*phi_s)*Diff.Dy2)...
                  + 2*m/Cn*diag(3*phi_s.^2-1)*Diff.Lap ...
                  - Cn*m*Diff.Lap2;
        
        [g,Jg] = get_g(phi_s);
        
        % Boundary condition at top and bottom boundaries
        topBottom       = Ind.bottom | Ind.top;
        y(topBottom)   = a_direction(topBottom,:)*(Diff.grad*phi_s) - g(topBottom);
        J(topBottom,:) = a_direction(topBottom,:)*Diff.grad-Jg(topBottom,:);
        
        % Boundary condition at left and right boundaries
        EYE = eye(M);
        y(Ind.left)    = phi_s(Ind.left) + 1;
        J(Ind.left,:)  = EYE(Ind.left,:);        
        
        y(Ind.right)    = phi_s(Ind.right) -1;
        J(Ind.right,:)  = EYE(Ind.right,:);
        
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