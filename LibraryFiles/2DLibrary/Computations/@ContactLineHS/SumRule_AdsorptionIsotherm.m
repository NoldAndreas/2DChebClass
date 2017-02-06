function [theta_sumrule,errTheta,I,errI] = SumRule_AdsorptionIsotherm(this,ST_LG) 

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
    
    ftMax = 13; 
    mark  = (this.AdsorptionIsotherm.FT < ftMax);
    
    mu_I = this.AdsorptionIsotherm.mu(mark) - mu_sat;
    ft_I = this.AdsorptionIsotherm.FT(mark);
            
    I = mu_I(1)*(ft_I(2)-ft_I(1))/2;
    for i = 2:(length(ft_I)-1)
        I = I + mu_I(i)*(ft_I(i+1)-ft_I(i-1))/2;
    end
    I = I + mu_I(end)*(ft_I(end)-ft_I(end - 1))/2;
    
    I = I + this.AdsorptionIsotherm.Coeff/(2*ftMax^2);
    
    I = I*(rhoLiq_sat-rhoGas_sat)/ST_LG;
    %(2) Compute contact angle from sum rule and compare with Young CA
    theta_sumrule = acos(1 - I)*180/pi;
    
    if(min(this.AdsorptionIsotherm.FT) < 0)
        theta_sumrule = 180 - theta_sumrule;
    end
    
    %fprintf(['Error for Contact Angle from Sum Rule for Disjoining Potential: ',num2str((theta_sumrule-this.alpha_YCA)*180/pi,2),' [deg]\n']);
    errI     = I-(1-abs(cos(this.alpha_YCA)));
    errTheta = errI/sin(theta_sumrule*pi/180)*180/pi;
    disp(['Integral expression: ',num2str(I),' +/- ',num2str(errI)]);
    fprintf(['Contact Angle from Sum Rule for Disjoining Potential: ',num2str(theta_sumrule),' [deg] +/_ ',num2str(errTheta),' [deg]\n' ]);
    
end