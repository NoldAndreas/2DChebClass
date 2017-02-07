function [contactlinePos,contactAngle_deg,y1Interface,ak]= ExtrapolateInterface(this,rho,y2,ak_ig,y2Lim,opts,t)
    if(nargin < 6)
        opts = {};    
    end
            
	%Based on values of the interface away from the wall at (y_1^k,y_2^k)
    %Interpolate to 
    % (EQ) y_1Interface = c*y_2+sum_k  a_k exp(-lambda_k *y_2^k)
    
    y2B     = 4;
    c       = 1/tan(this.optsNum.PhysArea.alpha_deg*pi/180);
    nk      = 1; %1;
    IP_1D   = barychebevalMatrix(this.IDC.Pts.x1,0);        
    rho0    = IP_1D*this.rho1D_lg;
    
    if(~IsOption(opts,'analyse'))             
    
       % y2k     = [4;5];
        %y2k     = (4:0.1:5.5)';    
        y2k_full     = (y2Lim(1):0.1:y2Lim(2))';        
        y1k_full     = zeros(size(y2k_full));
        
        lambdak = [0.5;0.6]; %Doesnt have to be equal y1k, I think.
        CartPts = this.IDC.GetCartPts;

        %1st step: get y1k                    
        fsolveOpts  = optimset('Display','off');        
        
        for k = 1:length(y2k_full)
            pt.y2_kv    = y2k_full(k);  
            
            [~,mark1] = min(abs(CartPts.y2_kv-pt.y2_kv));
            mark1     = (CartPts.y2_kv == CartPts.y2_kv(mark1));
            
            [~,mark2] = min(abs(rho(mark1)-rho0));
            y1I       = CartPts.y1_kv(mark1);
            y1I       = y1I(mark2);
                                                 
            y1k_full(k) = fsolve(@rhoX1,y1I,fsolveOpts);                                
        end    
            
        %2nd step: solve (EQ) for a_k
        if((nargin < 4) || isempty(ak_ig))
             if(nk == 2)
                ak_ig = [0;0.5;0;0.5];
             else
                ak_ig = [0;0.5];
             end         
        end
        
        y1k = [y1k_full(1),y1k_full(end)];
        y2k = [y2k_full(1),y2k_full(end)];
        ak_if  = fsolve(@f_error,ak_ig,fsolveOpts);                
        
        if(length(y2k) > 2)
            y1k = y1k_full;
            y2k = y2k_full;
            %'Display','final'
            %ak_ig = [ak_ig(1);sqrt(ak_ig(2))];
            %opts_min  = optimset('TolX',1e-12,'TolFun',1e-12,'MaxFunEvals',10000);
            opts_min  = optimset('TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',10000);
            %ak_ig     = [(y1k(end)-y1k(1))/(y2k(end)-y2k(1));mean(y1k)];
            ak        = fminsearch(@f_error_Av,ak_ig,opts_min);                
            %ak        = fminsearch(@f_error_Av_log,ak,opts_min);        
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
        y1Interface         = f_LinearExtrapolation(y2,ak);                 
    else  
        ak          = ak_ig;
        y1Interface = f_LinearExtrapolation(y2,ak);                 
        
        if(IsOption(opts,'subplot'))
            subplot(1,2,1);
            PlotDensityContours(this,rho);  hold on;      
            plot(y1Interface,y2,'linewidth',3);
        
            subplot(1,2,2);
            plot(y2,y1Interface - c*y2);
        else
            %**Get Isoline
            y2Plot      = (2:0.1:15)';
            y2Plot_Fit  = (0:0.1:7)';
            y1Plot      = zeros(size(y2Plot));
            fsolveOpts  = optimset('Display','off');
            for k = 1:length(y2Plot)
                pt.y2_kv    = y2Plot(k); 
                if(k>1)
                    y1Guess = y1Plot(k-1);
                else
                    y1Guess = 0;
                end
                    y1Plot(k)   = fsolve(@rhoX1,y1Guess,fsolveOpts);                                
            end    
            y1Plot_Fit = f_LinearExtrapolation(y2Plot_Fit,ak);
            
            if(~IsOption(opts,'full'))
                y1Plot     = y1Plot - c*y2Plot;
                y1Plot_Fit = y1Plot_Fit - c*y2Plot_Fit;
            end                
            
            plot(y1Plot,y2Plot,'k'); hold on;        
            if(IsOption(opts,'fit'))
                plot(y1Plot_Fit,y2Plot_Fit,'k:'); hold on;    
            end
            
            if(IsOption(opts,'full'))
                xlabel('$y_1$','Interpreter','Latex');
            else
                xlabel('$y_1 - y_2/\tan \theta$','Interpreter','Latex');
            end
            ylabel('$y_2$','Interpreter','Latex');      
            if(nargin > 6)
                text(y1Plot(end),y2Plot(end),...
                     ['$t = ',num2str(t,3),'$'],'Interpreter','Latex','rotation',90);     
            end
            

            %Compute R-Value in Interval [4,10]:
            %ybar = mean(y1Plot);
            %R = 1 - sum((y1Plot-y1Plot_Fit).^2)/sum((y1Plot-ybar).^2);
            %disp(['for t = ',num2str(t,3),', R = ',num2str(R)]);
        end
    end
    
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
        z_abs = sum((f_error(ak)).^2);
    end

    function z = f_error(ak)
        z      = f(y2k,ak) - y1k;
    end

    function z = f_error_Log(ak)
        z      = log(abs(f(y2k,ak) - y1k));
    end

%    function z = f_errorSq(ak)
%        z      = f_error([ak(1);ak(2).^2]);
%    end

    function [z,dz] = f(y2,ak)        
        z  = ak(1)*y2 + ak(2);
        dz = ak(1);
        %z      = c*y2  + exp(-y2*ak(2)).*(ak(1));
        %dz     = c     - ak(2)*ak(1)*exp(-y2*ak(2));        
        
        %if(nk == 2)
%            z      = z  + exp(-y2.^2*ak(4)).*(ak(3));
%            dz     = dz - 2*ak(3)*ak(4)*exp(-y2.^2*ak(4));        
%        end
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