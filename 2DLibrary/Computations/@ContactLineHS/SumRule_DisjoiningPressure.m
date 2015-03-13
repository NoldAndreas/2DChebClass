function errRel = SumRule_DisjoiningPressure(this,II_or_IV)

    if((nargin == 1) || strcmp(II_or_IV,'II'))
        anaDP   = this.disjoiningPressure_II;
    else
        anaDP   = this.disjoiningPressure_IV;
    end

    ST_LG = this.ST_1D.om_LiqGas;


    %******************************************
    %Check sum rule Eq. (11) in Henderson (2005)    
    %Integrate along y1, and check sum rule
    % - Int(DisjoiningPressure(y),y=-inf..inf) = \gammaLV*sin(theta)            

    Int_y1  = this.y1_SpectralLine.Int;
    sinGamm = sin(this.alpha_YCA)*ST_LG;
    err       = (Int_y1*anaDP + sinGamm);
    errRel    = err/sinGamm;            
    disp(['Integral expression: ',num2str(Int_y1*anaDP),' +/- ',num2str(err) , ' or relative: ',num2str(errRel*100),' percent']);

    h              = (-Int_y1*anaDP)/ST_LG;
    estTheta       = asin(h)*180/pi;
    if(IsDrying(this))
        estTheta = 180 - estTheta;
    end
    error_estTheta = 180/pi*sum(Int_y1)*max(abs(anaDP(1)),abs(anaDP(end)))/ST_LG/cos(estTheta*pi/180);%*(1/sqrt(1+h^2));
    disp(['Theta from Sum rule = ',num2str(estTheta),' [deg] +/- ',num2str(err/ST_LG/cos(estTheta*pi/180)*180/pi),' [deg]']);                            
    disp(['Error from numerical estimate : ',num2str(error_estTheta),' [deg]']);
    %disp(['Error from difference to integral estimate : ',num2str(err/ST_LG/cos(estTheta*pi/180)*180/pi),' [deg]']);

    %PrintErrorPos(180/pi*(estTheta-this.alpha_YCA),'Estimated contact angle through sum rule integrating disjoining pressure [percent]');
end  