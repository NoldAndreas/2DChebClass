function [u,Psi,p] =  GetHuhScriven_Solution(Pts,theta)

    r    = Pts.y1_kv;
    phi  = Pts.y2_kv;      
    
    if(length(theta)==1)
        theta = theta*ones(size(r));
    end

    %***********************
    % Coefficients
    %***********************
    
    C     = (theta.*(pi-theta)-(sin(theta)).^2+0.5*(2*theta-pi).*sin(2*theta))*pi;
    a1    = (pi*sin(theta)+pi*cos(theta).*(pi-theta))./C;
    a2    = (-0.5*pi*sin(2*theta) - pi*(pi-theta))./C; 
    b1    = (pi*(sin(theta)-theta.*cos(theta)))./C;
    b2    = (pi*(theta-0.5*sin(2*theta)))./C;

    %***********************
    % Functions 
    %***********************
    
    Psi_A = r.*(a1.*phi.*sin(phi-theta) + a2.*(phi-theta).*sin(phi));
    Psi_B = r.*(b1.*(phi-pi).*sin(phi-theta) + b2.*(phi-theta).*sin(phi));
    
    u_rA  = (a1.*(sin(phi-theta) + phi.*cos(phi-theta)) +...
                  a2.*(sin(phi) +(phi-theta).*cos(phi)));
    u_rB  = (b1.*(sin(phi-theta) + (phi-pi).*cos(phi-theta)) +...
                  b2.*(sin(phi) +(phi-theta).*cos(phi)));
              
    u_tA = -(a1.*phi.*sin(phi-theta) + a2.*(phi-theta).*sin(phi));
    u_tB = -(b1.*(phi-pi).*sin(phi-theta) + b2.*(phi-theta).*sin(phi));

    u_A        = [u_rA ; u_tA];
    u_B        = [u_rB ; u_tB];
    
    pA         = 2*(a1.*sin(phi-theta) + a2.*sin(phi))./r;
    pB         = 2*(b1.*sin(phi-theta) + b2.*sin(phi))./r;            

    %***********************
	% Merge
    %***********************
    
    mask_A     = Pts.y2_kv <= theta;
    mask_B     = Pts.y2_kv > theta;
    
    Psi         = zeros(size(r));
    Psi(mask_A) = Psi_A(mask_A);
    Psi(mask_B) = Psi_B(mask_B);

    mask_A2    = repmat(mask_A,2,1);
    mask_B2    = repmat(mask_B,2,1);
    
    u          = zeros(size([r;r]));
    u(mask_A2) = u_A(mask_A2);
    u(mask_B2) = u_B(mask_B2);
        
    p          = zeros(size(r));
    p(mask_A)  = pA(mask_A);
    p(mask_B)  = pB(mask_B);    

end