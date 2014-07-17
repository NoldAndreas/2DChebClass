function [rho,mu] = GetPointAdsorptionIsotherm(this,ell)

    %planarDisjoiningPressure = -(this.AdsorptionIsotherm_Mu-mu_sat)*(rhoLiq_sat-rhoGas_sat);

    ellAD = ell;%min(this.AdsorptionIsotherm_FT);
	IP    = barychebevalMatrix(this.AdsorptionIsotherm_FT,ellAD);
    
	mu    = IP*this.AdsorptionIsotherm_Mu;
    rho   = IP*this.AdsorptionIsotherm_rho;
    
	mu(ellAD<min(this.AdsorptionIsotherm_FT)) = 0;
	mu(ellAD>max(this.AdsorptionIsotherm_FT)) = 0;        
    
end