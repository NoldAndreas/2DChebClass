function theta = ComputeContactAngle(om_wallGas,om_wallLiq,om_LiqGas)
%Returns contact angle in rad.
%
%                      /
%                     /     Liquid
%      Gas           /  
%                   /theta
% ----------------------------------------
%               Wall
%
    ratio = (om_wallGas-om_wallLiq)/om_LiqGas;
    if(ratio > 1)
        theta = 0;
    elseif(ratio < -1)
        theta = pi;
    else
        theta = acos(ratio);
    end
    fprintf(['Contact angle from surface tensions = ',num2str(theta*180/pi),'[deg] \n']);
end