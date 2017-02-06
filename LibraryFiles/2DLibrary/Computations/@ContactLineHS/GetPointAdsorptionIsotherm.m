function [rho,mu] = GetPointAdsorptionIsotherm(this,ell)

    %planarDisjoiningPressure = -(this.AdsorptionIsotherm_Mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);
    ft    = abs(this.AdsorptionIsotherm.FT);

    ellAD = abs(ell);%min(this.AdsorptionIsotherm_FT);
	IP    = barychebevalMatrix(ft,ellAD);
    
	mu    = IP*this.AdsorptionIsotherm.mu;
    rho   = IP*this.AdsorptionIsotherm.rho;
    
	mu(ellAD<min(ft)) = 0;
	mu(ellAD>max(ft)) = 0;        
    
end