function [u] = GetSeppecherSolution(Pts,Vc,D_A,D_B,theta)

    r   = Pts.y1_kv;
    th  = Pts.y2_kv;      

    %1) Contribution from mass flux through the contact line
    u_rA1 = (1./r).*((D_A + D_B)/pi*(-1 + sin(theta-2*th)/sin(theta) + ...
                                 theta*sin(2*th)/(sin(theta)^2) )+...
                     - D_A*sin(2*th)/(sin(theta))^2); 
    u_rA1(u_rA1==inf)   = 0;
    u_rA1(u_rA1==-inf)  = 0;
    u_rA1(isnan(u_rA1)) = 0;
    u_rB1               = u_rA1;

    u_tA1      = zeros(size(u_rA1));
    u_tB1      = zeros(size(u_rA1));

    u_A        = [u_rA1 ; u_tA1];
    u_B        = [u_rB1 ; u_tB1];

    mask_A     = Pts.y2_kv <= theta;
    mask_B     = Pts.y2_kv > theta;

    mask_A2    = repmat(mask_A,2,1);
    mask_B2    = repmat(mask_B,2,1);

    u          = zeros(size(u_A));
    u(mask_A2) = u_A(mask_A2);
    u(mask_B2) = u_B(mask_B2);
    
    u          = u + Vc*GetHuhScriven_Solution(Pts,theta);   

end