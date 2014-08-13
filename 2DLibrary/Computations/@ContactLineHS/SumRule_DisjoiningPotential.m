function SumRule_DisjoiningPotential(this,ST_LG)

    if(nargin < 2)
        ST_LG = this.ST_1D.om_LiqGas;
    end
       
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