function LubricationDisjoiningPressure()

    %% Initialization
    shape    = struct('N',100,'L',10);
    PlotArea = struct('yMin',-20,'yMax',20,'N',200);    
    ST       = 1; %surfaceTension
    h0       = 1;
    Ca       = 0.01;
   % Theta    = 30*pi/180; % Angle of inclindation of substrate
%    Bond     = 10^-1; % Bond number

    hShape = InfSpectralLine(shape);
    hShape.ComputeAll(PlotArea);
    IntM   = hShape.IntegrationY();    
    y      = hShape.Pts.y;
    Diff   = hShape.Diff;
    
    [~,CA] = DisjoiningPressure(0,[]);
    disp(['Equilibrium contact angle: ',num2str(CA*180/pi),' [deg]'])
    hPInf  = tan(CA);
    
    %% Equilibrium
    %
    %
    % $$p = - \frac{d^2 h}{d x^2} - \frac{\Pi(h)}{ \gamma} = 0$$
    % 
    hP_iguess = hPInf*(1+tanh(y/2))/2;        
    hPEq  = fsolve(@f,hP_iguess);
    
    figure('name','Equilibrium');
    subplot(1,2,1);
    hShape.plot(hPEq);
    subplot(1,2,2);
    plot(y,h0+IntM*hPEq); xlim([-10 10]);
    %hShape.plot(h0+IntM*hP);
    
    %% Dynamics
    %
    % $$ Ca \frac{\partial h}{\partial t} =  \frac{1}{3} \frac{\partial}{\partial x}(h^3 \frac{\partial}{\partial x} p ) $$
    % 
    mM      = ones(size(hShape.Pts.y));
    mM(1)   = 0;    
    mM(end) = 0;
	opts    = odeset('RelTol',10^-10,'AbsTol',10^-10,'MvPattern',diag(mM));
    [outTimes,hP_t] = ode15s(@dhdt,0:0.2:10,hPEq,opts);
    hP_t = hP_t.';
    
    for i=1:length(outTimes)        
        hP  = hP_t(:,i);                
        h_  = h0 + IntM*hP;
        
        subplot(2,1,1);
        plot(y,h_); xlim([-30 10]);
        subplot(2,1,2);
        hShape.plot(hP);  
        title(['t = ', num2str(outTimes(i))]);
        zlim([0 max(hPEq)]);
        hold off;     
        pause(0.1);
        %Record(i,giffile,0.1);       
    end        
    
    
    
    %% Auxiliary Functions    
    function g = dhdt(t,hP)
        h     = h0+IntM*hP;
        hPP   = Diff.Dy*hP;
        hPPP  = Diff.DDy*hP;
        hPPPP = Diff.DDDy*hP;
        
        [dp,CA,ddp] = DisjoiningPressure(t,h);
        
        g = (h.^2).*(hPPP + ddp.*hP/ST) + (h.^3)/3.*(hPPPP + dp.*hPP/ST + ddp.*(hP.^2)/ST);
        g = -Diff.Dy*g/Ca;
        
        %Boundary Conditions:
        g(1)   = hP(1);
        g(end) = hP(end) - hPEq(end);
    end
    
    function y = f(hP)
        y        = staticCL(hP);
        y(1)     = hP(1);
      %  y(end)   = hP(end)-tan(CA);
        y(end) = IntM(end/2,:)*(hP-hP_iguess);%hP(end/2) - hP_iguess(end/2);
      %  y(end) = hP(end)-hP_iguess(end);
    end    
    
    function [dp,CA,ddp] = DisjoiningPressure(t,ell)        
       a   = 0.1*exp(-t);
       dp  = a*(1./ell.^9 - 1./ell.^3);
       ddp = a*(-9./ell.^10 + 3./ell.^4);
       CA = acos(1 - (3/8)*a/ST);       
    end

    function g = staticCL(hP)
        h = h0+IntM*hP;
        g = DisjoiningPressure(0,h) + ST*Diff.Dy*hP;% + ...
                %Bond*(h*cos(Theta)-y*sin(Theta));
    end

end

