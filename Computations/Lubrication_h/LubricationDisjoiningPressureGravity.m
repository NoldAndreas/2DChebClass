function LubricationDisjoiningPressureGravity()

    %% Initialization
    shape    = struct('N',200,'L',7);
    PlotArea = struct('yMin',-10,'yMax',10,'N',200);    
    ST       = 1; %surfaceTension
    h0       = 1;
    Theta    = 20*pi/180; % Angle of inclindation of substrate
    Bond     = 10^-2; % Bond number

    hShape = InfSpectralLine(shape);
    hShape.ComputeAll(PlotArea);
    IntM   = hShape.IntegrationY();    
    y      = hShape.Pts.y;
    
    [~,CA] = DisjoiningPressure([]);
    disp(['Equilibrium contact angle: ',num2str(CA*180/pi),' [deg]'])
    hPInf  = tan(CA);
    
    %% Computation
    hP_iguess = tan(Theta)*(1+tanh(y/2))/2; 
    
    opts = optimset('MaxFunEvals',10^6,'MaxIter',1000);
    hP   = fsolve(@f,hP_iguess);%,opts);
    
    subplot(1,2,1);
    hShape.plot(hP);
    subplot(1,2,2);
    plot(y,h0+IntM*hP); xlim([-10 10]);
    %hShape.plot(h0+IntM*hP);

    %% Auxiliary Functions
    function y = f(hP)
        y        = staticCL(hP);
        y(1)     = hP(1);
      %  y(end)   = hP(end)-tan(CA);
      %  y(end) = IntM(end/2,:)*(hP-hP_iguess);%hP(end/2) - hP_iguess(end/2);
        y(end) = hP(end)-hP_iguess(end);
    end    
    
    function [dp,CA] = DisjoiningPressure(ell)        
       a  = 0.1;
       dp = a*(1./ell.^9 - 1./ell.^3);
       CA = acos(1 - (3/8)*a/ST);       
    end

    function g = staticCL(hP)
        h =h0+ IntM*hP;
        g = DisjoiningPressure(h) + hShape.Diff.Dy*hP + ...
                Bond*((h-h0)*cos(Theta)-(1+tanh(y))/2*sin(Theta));
    end

end