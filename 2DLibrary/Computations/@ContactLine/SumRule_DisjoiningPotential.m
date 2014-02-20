function SumRule_DisjoiningPotential(this,ST_LG)

    if(nargin < 2)
        ST_LG = this.ST_1D.om_LiqGas;
    end
    
    %******************************************
    %Check sum rule Eq. (12) in Henderson (2005)
    %-int(DisjoiningPot,h=h_0..inf) = gamma_LV*(1-cos(theta))
    %(1) Integrate Disjoining Pressure    
    
    mu_sat         = this.optsPhys.mu_sat;    
    rhoLiq_sat     = this.optsPhys.rhoLiq_sat;
    rhoGas_sat     = this.optsPhys.rhoGas_sat;
    
    ftMax = 20;
    mark  = (this.AdsorptionIsotherm_FT < ftMax);
    
    mu_I = this.AdsorptionIsotherm_Mu(mark) - mu_sat;
    ft_I = this.AdsorptionIsotherm_FT(mark);
            
    I = mu_I(1)*(ft_I(2)-ft_I(1))/2;
    for i = 2:(length(ft_I)-1)
        I = I + mu_I(i)*(ft_I(i+1)-ft_I(i-1))/2;
    end
    I = I + mu_I(end)*(ft_I(end)-ft_I(end - 1))/2;
    
    I = I + this.AdsorptionIsotherm_Coeff/(2*ftMax^2);
    
    I = I*(rhoLiq_sat-rhoGas_sat);
    %(2) Compute contact angle from sum rule and compare with Young CA
    theta_sumrule = acos(1 - I/ST_LG);
    %fprintf(['Error for Contact Angle from Sum Rule for Disjoining Potential: ',num2str((theta_sumrule-this.alpha_YCA)*180/pi,2),' [deg]\n']);
    fprintf(['Contact Angle from Sum Rule for Disjoining Potential: ',num2str((theta_sumrule)*180/pi),' [deg]\n']);
    
    %******************************************
    %Check sum rule Eq. (11) in Henderson (2005)    
    %Integrate along y1, and check sum rule
    % - Int(DisjoiningPressure(y),y=-inf..inf) = \gammaLV*sin(theta)
    y1      = this.y1;
    Int_y1  = this.Int_y1;
    
    sinGamm = sin(this.alpha_YCA)*ST_LG;
    err     = (Int_y1*this.disjoiningPressure + sinGamm)/sinGamm*100;
    PrintErrorPos(err,'Normal force balance, error [percent]');
    
    h              = (-Int_y1*this.disjoiningPressure)/ST_LG;
    estTheta       = asin(h)*180/pi;
    error_estTheta = 180/pi*sum(Int_y1)*max(abs(this.disjoiningPressure(1)),abs(this.disjoiningPressure(end)))/ST_LG*(1/sqrt(1+h^2));
    disp(['Theta from Sum rule = ',num2str(estTheta),' [deg] /+- ',num2str(error_estTheta),' [deg]']);    
    
    %PrintErrorPos(180/pi*(estTheta-this.alpha_YCA),'Estimated contact angle through sum rule integrating disjoining pressure [percent]');

end