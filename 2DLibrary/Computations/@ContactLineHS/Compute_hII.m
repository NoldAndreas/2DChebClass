function Compute_hII(this)

    y1   = this.y1_SpectralLine.Pts.y;
    DP   = this.disjoiningPressure_II;
    ST   = this.ST_1D.om_LiqGas;
    
    if(this.optsNum.PhysArea.alpha_deg > 90)
        drying = true;
        DP     = flip(DP,1);
    else 
        drying = false;
    end
    
    %solves Eq. (17) from Henderson (2010)    
    %sets ell_2DDisjoiningPressure
    %uses disjoiningPressure,y1,Int_y1
    
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
    
    this.hII   = h_2D;
    
	function y = f3(h)
        y    = Diff.Dx*h-hP;
        y(1) = h(1);
    end

end