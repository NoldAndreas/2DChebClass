function Compute_hII(this,II_or_IV)

    y1   = this.y1_SpectralLine.Pts.y;
    if((nargin == 1) || strcmp(II_or_IV,'II'))
        DP   = this.disjoiningPressure_II;
    else
        DP   = this.disjoiningPressure_IV;
    end
        
    ST   = this.ST_1D.om_LiqGas;
    
    if(this.optsNum.PhysArea.alpha_deg > 90)
        drying = true;
        DP     = flip(DP,1);
    else 
        drying = false;
    end
    
    %solves Eq. (17) from Henderson (2010), (29) from Nold et.al. (2014)
    % 
    % h_II'(x) = -tan( asin ( 1/gamma \int[ Pi_II(x'),x'=-infinity..x] ))
    % 
    
    
    Diff = barychebdiff(y1,2);    
    
    %**Alternative Computation
	IntMY1 = zeros(length(y1));
	vh     = zeros(1,length(y1));
    for i = 2:length(y1)
        %Matrix for integral int(v(y),y=y1..Y)
        h           = y1(i) - y1(i-1);
        vh([i-1,i]) = vh([i-1,i]) + h/2;
        IntMY1(i,:)   = vh;
    end
    
    hP   = -tan(asin(IntMY1*DP/ST));
    h_2D = fsolve(@f3,zeros(size(y1)));      
    if(drying)
        h_2D     = flip(h_2D,1);
    end
    
    if((nargin == 1) || strcmp(II_or_IV,'II'))
        this.hII   = h_2D;
    else
        this.hIV   = h_2D;
    end
        
    
	function y = f3(h)
        y    = Diff.Dx*h-hP;
        y(1) = h(1);
    end

end