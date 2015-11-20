function [contactlinePos,contactAngle_deg,y1Interface,ak]= ExtrapolateInterface(this,rho,y2,ak_ig,y2Lim)
            
	%Based on values of the interface away from the wall at (y_1^k,y_2^k)
    %Interpolate to 
    % (EQ) y_1Interface = c*y_2+sum_k  a_k exp(-lambda_k *y_2^k)
    
    
   % y2k     = [4;5];
    %y2k     = (4:0.1:5.5)';    
    y2k     = (y2Lim(1):0.1:y2Lim(2))';
    nk      = 1; %1;
    
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
         if(nk == 2)
            ak_ig = [0;0.5;0;0.5];
         else
            ak_ig = [0;0.5];
         end
         
%         if(c > 1)
%             ig = [1;0.5];
%         else
%             ig = [-1;0.5];
%         end
    end
        
%    if(ak_ig(2) < 0)
%        ak_ig = [1;0.5];
    %end
    
    if(length(y2k) > 2)
        %'Display','final'
        %ak_ig = [ak_ig(1);sqrt(ak_ig(2))];
        opts  = optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',10000);
        ak    = fminsearch(@f_error_Av,ak_ig,opts);
        %ak    = [ak(1);ak(2).^2];
        if(ak(2) < 0)
            ak(1) = fminsearch(@f_error_Av_Zero,0,opts);
            ak(2) = 0;
        end
    else
        ak  = fsolve(@f_error,ak_ig,fsolveOpts);
    end
    
	disp(['Fitting coefficient ak = ',num2str(ak'),' with |error| = ',num2str(f_error_Av(ak))]);
    

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

    function z_abs = f_error_Av_Zero(ak)
        z_abs = sum(abs(f_error([ak;0])));
    end

    function z_abs = f_error_Av(ak)
        z_abs = sum(abs(f_error(ak)));
    end

    function z = f_error(ak)
        z      = f(y2k,ak) - y1k;
    end

%    function z = f_errorSq(ak)
%        z      = f_error([ak(1);ak(2).^2]);
%    end

    function [z,dz] = f(y2,ak)        
        z      = c*y2  + exp(-y2*ak(2)).*(ak(1));
        dz     = c     - ak(2)*ak(1)*exp(-y2*ak(2));        
        
        if(nk == 2)
            z      = z  + exp(-y2.^2*ak(4)).*(ak(3));
            dz     = dz - 2*ak(3)*ak(4)*exp(-y2.^2*ak(4));        
        end
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