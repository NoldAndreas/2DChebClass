function [u] = GetSeppecherSolution(Pts,Vc,D_A,D_B,Phi) %,p,mask_A,mask_B]

    r   = Pts.y1_kv;
    th  = Pts.y2_kv;      

    %1) Contribution from mass flux through the contact line
    u_rA1 = (1./r).*((D_A + D_B)/pi*(-1 + sin(Phi-2*th)/sin(Phi) + ...
                                 Phi*sin(2*th)/(sin(Phi)^2) )+...
                     - D_A*sin(2*th)/(sin(Phi))^2); 
    u_rA1(u_rA1==inf)   = 0;
    u_rA1(u_rA1==-inf)  = 0;
    u_rA1(isnan(u_rA1)) = 0;
    u_rB1               = u_rA1;

    u_tA1 = 0;
    u_tB1 = 0;

    %2) Contribution from the Huh and Scriven solution
    C     = (Phi*(pi-Phi)-(sin(Phi))^2+0.5*(2*Phi-pi)*sin(2*Phi))    *pi;
    a1    = pi*sin(Phi)+pi*cos(Phi)*(pi-Phi);
    a2    = -0.5*pi*sin(2*Phi) - pi*(pi-Phi); 
    b1    = pi*(sin(Phi)-Phi*cos(Phi));
    b2    = pi*(Phi-0.5*sin(2*Phi));

    u_rA2 = Vc/C*(a1*(sin(th-Phi) + th.*cos(th-Phi)) +...
                  a2*(sin(th) +(th-Phi).*cos(th)));
    u_rB2 = Vc/C*(b1*(sin(th-Phi) + (th-pi).*cos(th-Phi)) +...
                  b2*(sin(th) +(th-Phi).*cos(th)));

    u_tA2 = -Vc/C*(a1*th.*sin(th-Phi) + a2*(th-Phi).*sin(th));
    u_tB2 = -Vc/C*(b1*(th-pi).*sin(th-Phi) + b2*(th-Phi).*sin(th));

    u_rA  = u_rA1 + u_rA2;
    u_tA  = u_tA1 + u_tA2;

    u_rB  = u_rB1 + u_rB2;
    u_tB  = u_tB1 + u_tB2;             

    u_A   = [u_rA ; u_tA];
    u_B   = [u_rB ; u_tB];

    mask_A    = Pts.y2_kv < Phi;
    mask_B    = Pts.y2_kv > Phi;

    mask_A2    = repmat(mask_A,2,1);
    mask_B2    = repmat(mask_B,2,1);

    u          = zeros(size(u_A));
    u(mask_A2) = u_A(mask_A2);
    u(mask_B2) = u_B(mask_B2);

%    p  = zeros(size(r));
 %   if(((D_A~= 0)||(D_B~=0)))
%        disp('GetSeppecherSolution: Pressure for DA,DB ~= 0 not yet implemented.');
%    else
     %   pA = 2./r.*(a1*sin(th-Phi) + a2*sin(th));
%        pB = 2./r.*(b1*sin(th-Phi) + b2*sin(th));            
%        p(mask_A) = pA(mask_A);
%        p(mask_B) = pB(mask_B);
%    end            

end