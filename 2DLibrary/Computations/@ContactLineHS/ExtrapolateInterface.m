function [contactlinePos,contactAngle_deg,y1Interface]= ExtrapolateInterface(this,rho,y2)
            
	%Based on values of the interface away from the wall at (y_1^k,y_2^k)
    %Interpolate to 
    % (EQ) y_1Interface = c*y_2*(1+sum_k  a_k exp(-lambda_k *y_2^k))
       
    y2k     = [4;5];
    lambdak = [0.5;0.6]; %Doesnt have to be equal y1k, I think.
    c       = tan(this.optsNum.PhysArea.alpha_deg*pi/180);
    CartPts = this.IDC.GetCartPts;
        
    %1st step: get y1k    
    IP_1D   = barychebevalMatrix(this.IDC.Pts.x1,0);
    rho0    = IP_1D*this.rho1D_lg;
    y1k     = zeros(size(y2k));
    fsolveOpts  = optimset('Display','off');
    
    for k = 1:length(y2k)
        pt.y2_kv    = y2k(k);        
        y1k(k)      = fsolve(@rhoX1,0,fsolveOpts);                                
    end
    
    %2nd step: solve (EQ) for a_k
    ig     = ones(size(y2k));
    %ig = [1;1];
    ak     = fsolve(@f_error,ig);
    %ak

    %3rd step: set interface, contact line position
    y2min          = min(CartPts.y2_kv);
    [contactlinePos,dz] = f(3,ak);   % 1+y2min
    contactAngle_deg    = atan(dz)*180/pi;
    if(nargin > 2)
        y1Interface = f(y2,ak);
    end

    PlotDensityContours(this,rho);  hold on;      
    plot(y1Interface,y2,'linewidth',3)
     
    function z = rhoX1(y1)
        pt.y1_kv = y1;
        IP       = this.IDC.SubShapePtsCart(pt);
        z        = IP*rho-rho0;
    end

    function z = f_error(ak)
        z      = f(y2k,ak) - y1k;
    end

    function [z,dz] = f(y2,ak)
        %z      = c*y2 + exp(-y2*lambdak')*ak;
        %dz     = c - (exp(-y2*lambdak')*(ak.*lambdak));
        
        z      = c*y2 + exp(-y2*0.5).*(ak(1)+ak(2).*(log(1./(0.5+y2))).^2) ;
        %z      = c*y2 + exp(-y2*0.5).*(ak(1)+ak(2)*(1./(0.5+y2)));
        dz     = 0;
        %dz     = c - (exp(-y2*lambdak')*(ak.*lambdak));
        
        %lambdak = 2;
        %z      = c*y2 + (1./(1+y2).^(lambdak))*ak;        
        %dz     = c - lambdak*(1./(1+y2).^(lambdak+1))*ak;
        
        %lambdak = 1;
        %z      = c*y2 + (1./(1+y2*lambdak'))*ak;        
        %dz     = c - (1./(1+y2*lambdak'))*(ak.*lambdak);
    end
     
end