function LubricationForcedWetting_old()
%% Reproduce computation in Fig. 1 of Eggers PoF,17,082106 (2005)

%% ODE to solve
% 
% $$\frac{-\delta}{h^2+\lambda h} = h'''-h'+ \theta$$
% 
    %% Parameters
    delta   = 0.3;
    lambda3 = 1e-5;
    theta   = 1;
        
    N = 500;

    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape      = struct('N',500,'L',1);    
    plotShape  = struct('yMin',0,'yMax',4,'N',N);

    HIS        = HalfInfSpectralLine(shape);    
    [Pts,Diff,Int] = HIS.ComputeAll(plotShape);            
    y    = Pts.y; 
    DDy  = Diff.DDy;      
    
    IntM = zeros(length(y));
	vh     = zeros(1,length(y));
    for i = 2:length(y)
        %Matrix for integral int(v(y),y=y..Y)
        hh          = y(i) - y(i-1);
        vh([i-1,i]) = vh([i-1,i]) + hh/2;
        IntM(i,:)    = vh;
    end
    
    
	shapeLine  = struct('yMin',0,'yMax',1,'N',300);    
	SL         = SpectralLine(shapeLine);
    SL.ComputeAll();
	IntL       = SL.Int;
    
    %% Solve ODE
    hP  = fsolve(@ODE,ones(N,1));
    
    %% Get inner and outer region results
    % Inner Region    
    L         = lambda3/exp(1);%*(1-(pi^2-1)*delta/6);
    hIn3Asymp = 1 + 3*delta*log(y/L);
    % Outer Region    
    thetaAP = ComputeThetaAP();
    
    C          = ComputeC(theta,thetaAP);
    c          = exp(C);
%    thetaAP    = (1-3*delta*log(c*L))^(1/3);%2.32;%
    disp(['thetaAP = ',num2str(thetaAP)]);
    
    % h0        = theta*y+(theta-thetaAP)*(exp(-y)-1);
    h0P         = theta-(theta-thetaAP)*exp(-y);
    hOut3Asymp  = thetaAP^3 + 3*delta*log(y) + 3*delta*C;    
    
    % h1P(y): derivative of first order outer solution (see Eq. (16))
    
    % f(y) = -1/(h0(y))^2    
    %% 1st term
    % 
    % $$-int(exp(x-t)*f(t),t=x..infinty)/2$$
    % 
    T1 = zeros(N,1);
    for i = 1:N
        t     = y;
        T1(i) = -Int*(exp(-t).*f(y(i)+t))/2;
    end
    
    %% 2nd term
    % 
    % $$-\frac{1}{2}\int_1^x e^{t-x} f(t)dt$$
    % 
	T2 = zeros(N,1);
    for i = 1:(N-1)
        t     = 1 + SL.Pts.y*(y(i)-1);
        t(1)  = 1;
        T2(i) = (1-y(i))/2*IntL*(exp(t-y(i)).*f(t));
    end        
    T2(N) = 0;
        
    %% 3rd term: -K*exp(-x)   
    % 
    % $$K = \int_1^\infty (e^{-t}/2 -1) f(t) dt + \int_0^1 (\cosh(t)-1) f(t)dt$$
    %    
    % First term
    t  = y + 1;    
    K  = Int*((exp(-t)/2-1).*f(t));
    
    % 2nd term
    tt        = SL.Pts.y;    
    integrand = (tt.^2).*f(tt)/2;
	integrand(1) = -1/(2*thetaAP^2);
	K = K + IntL*integrand;
    
    for n = 4:2:100
        integrand = (tt.^n).*f(tt)/factorial(n);
        integrand(1) = 0;
        K = K + IntL*integrand;
    end         
    
    T3 = -K*exp(-y);
    
    h1P = T1 + T2 + T3;

    %% Plot results
        
    h = figure('color','white');
	shift = struct('xmin',0,'xmax',4,'ymin',0,'ymax',2.5,'yref',0,'xref',0);
	PlotBackgroundImage(['Computations' filesep 'PDEs' filesep 'EggersPoF2005_Fig3_cut.gif'],shift);    
    HIS.plot(hP);   
    %Outer solution:
    plot(y,thetaAP+(delta/thetaAP^2)*(log(y*c)),'r','linewidth',1.5); %C/thetaAP^2+
    plot(y,1 + delta*log(y/L),'g','linewidth',1.5);
    plot(y,(thetaAP^3+3*delta*log(y*c)).^(1/3),'b','linewidth',1.5);
    plot(y,(1+3*delta*log(y/L)).^(1/3),'m--','linewidth',1.5);
    plot(y,h0P+delta*h1P,'k-.','linewidth',1.5);
    plot(y,h0P,'k:','linewidth',1.5);
    
    %plot(y,h0P+delta*h1P,'g:','linewidth',1.5);
    ylim([0 2.5]);
    xlim([0 4]);
    
    %copyobj(gcf,0);
    set(gca,'XScale','log');

    %% Right hand side of ODE
    function y = ODE(hP)
        
        h = IntM*hP;        
        y = delta./(h.^2+lambda3*h) + DDy*hP - hP + theta;
        
        % Boundary conditions 
        y(1)   = hP(1)-1;
        y(end) = hP(end)-theta;
        
    end

    function [y,coeff0_0,coeff1_0] = f(x,th,thAP)
       if(nargin < 3)
           th   = theta;
           thAP = thetaAP;
       end
       hh0   = th*x+(th-thAP)*(exp(-x)-1);
       y     = -1./hh0.^2;
       
       % second coefficient for expansion around x = 0
       coeff0_0 = -1/thAP^2;
       coeff1_0 = (th-thAP)/thAP^3;
    end

    function thAP = ComputeThetaAP()
        thAP = fsolve(@eqn1,theta);
    end

    function y = eqn1(thAP)
        y = 1 - 3*delta*(log(L)+ComputeCout(theta,thAP))-thAP^3;
        end
        
    function C = ComputeC(th,th_ap)
        
        %analytical expression for th=th_ap = thetaAP
%        C         = exp(-1)-expint(1) -1;         
%        for nn = 2:100
%            C = C - (-1)^(nn-1)/(factorial(nn)*(nn-1));
%        end               
        
        ytilde     = 1 + y;
        integrand1 = (1-exp(-ytilde)).*f(ytilde,th,th_ap)*(th_ap^2);
        C          = Int*integrand1;       
        
        ytilde                 = SL.Pts.y;
        [fh,coeff0_0,coeff1_0] = f(ytilde,th,th_ap);
        integrand2             = (1-exp(-ytilde)).*fh*th_ap^2 + 1./ytilde;
        integrand2(1)          = (coeff1_0-coeff0_0/2)*th_ap^2;
        C                      = C + IntL*integrand2;

        
    end


     function C = ComputeCout(th,th_ap)
        
        %analytical expression for th=th_ap = thetaAP
%        C         = exp(-1)-expint(1) -1;         
%        for nn = 2:100
%            C = C - (-1)^(nn-1)/(factorial(nn)*(nn-1));
%        end               
        
        ytilde     = 1 + y;
        integrand1 = (1-exp(-ytilde)).*f(ytilde,th,th_ap);         
        
        ytilde                 = SL.Pts.y;
        [fh,coeff0_0,coeff1_0] = f(ytilde,th,th_ap);
        integrand2             = (1-exp(-ytilde)).*fh+ 1./((th_ap^2)*ytilde);
        integrand2(1)          = coeff1_0-coeff0_0/2;
                
        C   = (th_ap^2)*(Int*integrand1 + IntL*integrand2) +2;

    end

end

