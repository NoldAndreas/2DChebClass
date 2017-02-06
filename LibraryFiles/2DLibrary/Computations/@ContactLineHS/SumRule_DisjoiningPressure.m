function [errRel,estTheta,error_estTheta,I,err] = SumRule_DisjoiningPressure(this,II_or_IV)

    if((nargin == 1) || strcmp(II_or_IV,'II'))
        anaDP   = this.disjoiningPressure_II;
        fprintf('II:\n');
    else
        anaDP   = this.disjoiningPressure_IV;
        fprintf('I:\n');
    end

    ST_LG = this.ST_1D.om_LiqGas;


    %******************************************
    %Check sum rule Eq. (11) in Henderson (2005)    
    %Integrate along y1, and check sum rule
    % - Int(DisjoiningPressure(y),y=-inf..inf) = \gammaLV*sin(theta)            

    Int_y1  = this.y1_SpectralLine.Int;    
   
    I         = (-Int_y1*anaDP)/ST_LG;    
    err       = I - sin(this.alpha_YCA);
    errRel    = err/I;            
    disp(['Integral expression: ',num2str(I),' +/- ',num2str(err) , ' or relative: ',num2str(errRel*100),' percent']);

    
    estTheta       = asin(I)*180/pi;
    if(IsDrying(this))
        estTheta = 180 - estTheta;
    end
    
    error_estTheta = err/ST_LG/cos(estTheta*pi/180)*180/pi;
    disp(['Theta from Sum rule = ',num2str(estTheta),' [deg] +/- ',num2str(error_estTheta),' [deg]']);                            
    
    error_estTheta_num = 180/pi*sum(Int_y1)*max(abs(anaDP(1)),abs(anaDP(end)))/ST_LG/cos(estTheta*pi/180);%*(1/sqrt(1+h^2));
    disp(['Error from numerical estimate : ',num2str(error_estTheta_num),' [deg]']);
    %disp(['Error from difference to integral estimate : ',num2str(err/ST_LG/cos(estTheta*pi/180)*180/pi),' [deg]']);

    %PrintErrorPos(180/pi*(estTheta-this.alpha_YCA),'Estimated contact angle through sum rule integrating disjoining pressure [percent]');
end  