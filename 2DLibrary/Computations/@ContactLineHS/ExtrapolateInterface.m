function [contactlinePos,contactAngle_deg,y1Interface,ak]= ExtrapolateInterface(this,rho,y2,ak_ig)
            
	%Based on values of the interface away from the wall at (y_1^k,y_2^k)
    %Interpolate to 
    % (EQ) y_1Interface = c*y_2*(1+sum_k  a_k exp(-lambda_k *y_2^k))
       
    y2k     = [4;5];
    lambdak = [0.5;0.6]; %Doesnt have to be equal y1k, I think.
    y2B     = 4;
    c       = 1/tan(this.optsNum.PhysArea.alpha_deg*pi/180);
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
    if((nargin < 4) || isempty(ak_ig))
         ak_ig = [0;0.5];
%         if(c > 1)
%             ig = [1;0.5];
%         else
%             ig = [-1;0.5];
%         end
    end
    ak     = fsolve(@f_error,ak_ig,fsolveOpts);    
	disp(['Fitting coefficient ak = ',num2str(ak')]);
    

    %3rd step: set interface, contact line position    
    y2Min = min(this.IDC.GetCartPts.y2_kv);
    [contactlinePos,dz] = f_LinearExtrapolation(y2Min,ak);
    contactAngle_deg    = mod(atan(1/dz)*180/pi,180);
    if(nargin > 2)
        %y1Interface = f(y2,ak);
        y1Interface = f_LinearExtrapolation(y2,ak);         
    end

    PlotDensityContours(this,rho);  hold on;      
    plot(y1Interface,y2,'linewidth',3)
    
%     y2kAna = [4:0.2:15];
%     for k = 1:length(y2kAna)
%         pt.y2_kv    = y2kAna(k);        
%         y1kAna(k)   = fsolve(@rhoX1,0,fsolveOpts);                                
%     end   
    
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
        
        %z      = c*y2 + exp(-y2*0.5).*(ak(1));
        %dz     = c - 0.5*ak(1)*exp(-y2*0.5);
        
        z      = c*y2 + exp(-y2*ak(2)).*(ak(1));
        dz     = c - ak(2)*ak(1)*exp(-y2*ak(2));
        
        %cexpL  = 2;
        %z      = c*y2 + exp(-y2*0.5).*(ak(1)+ak(2).*(log(1./(0.5+y2))).^cexpL) ;        
        %dz     = c - 0.5*exp(-y2*0.5).*(ak(1)+ak(2).*(log(1./(0.5+y2))).^cexpL)...
%                      + exp(-y2*0.5).*(-cexpL*ak(2).*((-log(0.5+y2)).^(cexpL-1))./(0.5+y2));
        
        %z      = c*y2 + exp(-y2*0.5).*(ak(1)+ak(2)*(1./(0.5+y2)));
        %dz     = c - (exp(-y2*lambdak')*(ak.*lambdak));
        
        %lambdak = 2;
        %z      = c*y2 + (1./(1+y2).^(lambdak))*ak;        
        %dz     = c - lambdak*(1./(1+y2).^(lambdak+1))*ak;
        
        %lambdak = 1;
        %z      = c*y2 + (1./(1+y2*lambdak'))*ak;        
        %dz     = c - (1./(1+y2*lambdak'))*(ak.*lambdak);
    end

    function [z,dz] = f_LinearExtrapolation(y2,ak)
        z           = zeros(size(y2));
        dz          = zeros(size(y2));
        
        mark        = (y2 >= y2B);
        [z(mark),dz(mark)]  = f(y2(mark),ak);
        
        [y1B,dy1B] = f(y2B,ak);        
        
        z(~mark)    = dy1B *(y2(~mark)-y2B) + y1B;
        dz(~mark)   = dy1B;
    end
     
end