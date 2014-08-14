function Compute_hI(this)
    
    fT   = this.AdsorptionIsotherm.FT;
	IntM = zeros(length(fT));
	vh   = zeros(1,length(fT));
    for i = 2:length(fT)        
        %Matrix for integral int(v(y),y=y1..Y)
        h           = fT(i) - fT(i-1);
        vh([i-1,i]) = vh([i-1,i]) + h/2;
        IntM(i,:)   = vh;
    end
    if(drying)
        fT = -fT;
    end
       
    planarDisjoiningPressure = -(this.AdsorptionIsotherm.mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);
    IntPot = -IntM*planarDisjoiningPressure; 
    
   % plot(fT,planarDisjoiningPressure); hold on;
%    plot(fT,IntPot,'m');
    checkSumRule = this.ST_1D.om_LiqGas*(1-cos(this.alpha_YCA));
 %   plot([min(fT) max(fT)],[checkSumRule checkSumRule],'k--');
    
    x         = ClenCurtFlip(150);
    CompDiff  = barychebdiff(x,2);    
    
    if(this.optsNum.PhysArea.alpha_deg == 60)
        hmax =  this.filmThickness(75)-min(this.filmThickness);    
        L = this.y1(75)-min(this.y1);
        z = (x+1)/2*L -L+y1(75);
    else
        hmax =  this.filmThickness(end)-min(this.filmThickness);    
        L = max(this.y1)-min(this.y1);
        z = (x+1)/2*L -L+max(y1);
    end
    DPhys.Dx  = CompDiff.Dx*(2/L);    
    ST_LG = this.ST_1D.om_LiqGas;
    h  = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    %this.ell_1DDisjoiningPressure   = h; 
    
    
    hmax = h_2D(end);
    L2 = 10;
    z2 = (x+1)/2*L2 -L2 +max(this.y1);
    DPhys.Dx  = CompDiff.Dx*(2/L2);        
    h2   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    hmax = h2(1);
    L3 = 12;
    z3 = (x+1)/2*L3 - L3 +min(z2);
    DPhys.Dx  = CompDiff.Dx*(2/L3);
    h3   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    hmax = h3(1);
    L4 = 12;
    z4 = (x+1)/2*L4 - L4 +min(z3);
    DPhys.Dx  = CompDiff.Dx*(2/L4);
    h4   = fsolve(@f6,hmax*ones(size(z)));%exp(-z));%ell_h);
    
    this.hI   = [h;h2;h3;h4];
    this.y1_I = [z;z2;z3;z4];
    
    function y = f6(h)
        dh = ((1-interfacePotential(h)/ST_LG).^(-2)-1).^0.5;
        y  = DPhys.Dx*h - dh;
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