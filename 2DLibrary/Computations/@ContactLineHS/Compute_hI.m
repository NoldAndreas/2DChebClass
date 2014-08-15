function Compute_hI(this)
    
    fT             = abs(this.AdsorptionIsotherm.FT);
    y1             = this.y1_SpectralLine.Pts.y;
    ST_LG          = this.ST_1D.om_LiqGas;      
    Dx  = 0;    
    N   = 100;    
    L   = 3;
           
    [~,IntPot] = GetDisjoiningPressure_I(this);    
    
%     if(~isempty(this.hIII))
%         [hmax,iM] = max(this.hIII);
%         hmax      = hmax - min(this.hIII);
%         y10       = y1(iM);    
%     else
        hmax = 15; y10 = 0;
%    end
    
    yN = 0; hN = hmax;    
    yC = []; hC = [];
    
    while(max(abs(imag(hN))) == 0)
        yC = [yN;yC];
        hC = [hN;hC];
        [yN,hN] = Compute_hI_Interval([yC(1)-L,yC(1)],hC(1));
    end
        
    this.hI   = hC;   
    
    if(IsDrying(this))        
        this.y1_I = y10 - yC;
    else
        this.y1_I = y10 + yC;
    end
    
    %    checkSumRule = this.ST_1D.om_LiqGas*(1-cos(this.alpha_YCA));
    
    function [y1,h2] = Compute_hI_Interval(interval,hRight)
        hmax      = hRight;
        SL        = SpectralLine(struct('yMin',interval(1),...
                                        'yMax',interval(2),'N',N));
        SL.ComputeAll();
        y1        = SL.Pts.y;
        Dx        = SL.Diff.Dy;
        h2        = fsolve(@f6,hmax*ones(N,1));%exp(-z));%ell_h);
        
    end
    
    function y = f6(h)
        dh = ((1-interfacePotential(h)/ST_LG).^(-2)-1).^0.5;
        y  = Dx*h - dh;
        y(end) = h(end) - hmax;
    end

    function dP1D = interfacePotential(ell) 
        
        ellAD = min(fT)+ell;
        IP    = barychebevalMatrix(fT,ellAD);
        dP1D  = IP*IntPot;
        dP1D(ellAD<min(fT)) = 0;
     %   dP1D(ellAD>max(this.AdsorptionIsotherm_FT)) = 0;
    end

end