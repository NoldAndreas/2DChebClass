function LubricationForcedWetting()
%% Reproduce computation in Fig. 1 of Eggers PoF,17,082106 (2005)

%% ODE to solve
% 
% $$\frac{-\delta}{h^2+\lambda h} = h'''-h'+ \theta$$
% 
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'ForcedWetting'],'ORG');    
    close all;
    
    %% Parameters
    delta   = 0.3;
    lambda  = 1e-5;
    theta   = 1;           

    N = 500;
    
    %% Initialization
    %     
    % # for computation of similarity solution
    % # for 2D area for plotting of streamfunction
    shape      = struct('N',N,'L',1);    
    plotShape  = struct('yMin',0,'yMax',4,'N',N);

    HIS            = HalfInfSpectralLine(shape);    
    [Pts,Diff,Int] = HIS.ComputeAll(plotShape);                
    y              = Pts.y; 
    DDy            = Diff.DDy;      
    
    IntM = zeros(length(y));
	vh     = zeros(1,length(y));
    for i = 2:length(y)
        %Matrix for integral int(v(y),y=y..Y)
        hh          = y(i) - y(i-1);
        vh([i-1,i]) = vh([i-1,i]) + hh/2;
        IntM(i,:)    = vh;
    end
    
    
	shapeLine  = struct('yMin',0,'yMax',1,'N',N);
	SL         = SpectralLine(shapeLine);
    SL.ComputeAll();
	IntL       = SL.Int;
    
    %% Solve ODE
     hP  = fsolve(@ODE,ones(N,1));
%     
%     fig1 = figure('color','white');
% 	shift = struct('xmin',0,'xmax',4,'ymin',0,'ymax',2.5,'yref',0,'xref',0);
% 	PlotBackgroundImage(['Computations' filesep 'PDEs' filesep 'EggersPoF2005_Fig3_cut.gif'],shift);    
%     HIS.doPlots(hP);            
%     ylim([0 2.5]);
%     xlim([0 4]); 
%     xlabel('$x$','Interpreter','Latex');
%     ylabel('$h$','Interpreter','Latex');
%     
%     print2eps([dirData filesep 'Eggers1'],fig1);
%     saveas(fig1,[dirData filesep 'Eggers1.fig']);
    
    
    %% Inner region results
    %L         = lambda/exp(1);%*(1-(pi^2-1)*delta/6);    
    Cin     = 1;%-log(1-(pi^2-1)*delta/6);
    
    th = 1:0.01:2.5;
    for i = 1:length(th)
        hh(i) = eqn1(th(i));
    end
    plot(th,hh);
        
    
    % Outer Region    
    thetaAP     = ComputeThetaAP();        
    Cout        = ComputeCout(theta,thetaAP);    
%    thetaAP    = (1-3*delta*log(c*L))^(1/3);%2.32;%
    disp(['thetaAP = ',num2str(thetaAP)]);
        
    h0P         = theta-(theta-thetaAP)*exp(-y);                    
    st = DataStorage('NextOrder',@ComputeNextOrder,struct('f','f1'),struct('f',@f1),false);
    h1 = st.hN; h1P = st.hNP;
    st = DataStorage('NextOrder',@ComputeNextOrder,struct('f','f2'),struct('f',@f2),false);
    h2 = st.hN; h2P = st.hNP;
    %[h1,h1P]    = ComputeNextOrder(@f1);        
    %[h2,h2P]    = ComputeNextOrder(@f2);

    %% Plot results
        
    fig1 = figure('color','white','Position',[0 0 800 800]);
	HIS.doPlots(hP);    
    %yG = y(y>lambda/3);
    yP = 10.^((-5.6:0.01:1)');
    
   % plot(y,thetaAP*ones(size(y)),'r:','linewidth',1.5); 
    
    plot(y,h0P,'k:','linewidth',1.5);
    plot(y,h0P+delta*h1P,'k--','linewidth',1.5);    
    plot(y,h0P+delta*h1P+delta^2*h2P,'k-.','linewidth',1.5);
    
    xlabel('$x$','Interpreter','Latex','fontsize',20);
    ylabel('$dh/dx$','Interpreter','Latex','fontsize',20);
    
    print2eps([dirData filesep 'Eggers1'],fig1);
    saveas(fig1,[dirData filesep 'Eggers1.fig']);
    
    ylim([0 2.5]);
    xlim([0 4]);             

    plot(yP,thetaAP+(delta/thetaAP^2)*(log(yP)+Cout),'r:','linewidth',1.5); 
    plot(yP,(thetaAP^3+3*delta*(log(yP)+Cout)).^(1/3),'b','linewidth',1.5);    
    
    Tout_overlap = (thetaAP^3+3*delta*(log(yP))).^(1/3) + delta*Cout*(thetaAP^3+3*delta*(log(yP))).^(-2/3);    
    plot(yP,Tout_overlap,'b:','linewidth',1.5);
    
    plot(yP,1 + delta*(log(yP/lambda)+Cin),'r:','linewidth',1.5);    
    plot(yP,(1+3*delta*(log(yP/lambda)+Cin)).^(1/3),'m--','linewidth',1.5);                    
    
    
    Tin_overlap = (1+3*delta*(log(yP/lambda))).^(1/3) + delta*Cin*((1+3*delta*(log(yP/lambda))).^2).^(-1/3);    
    plot(yP,Tin_overlap,'m:','linewidth',1.5);
    
    
    %copyobj(gcf,0);
    set(gca,'XScale','log');

    print2eps([dirData filesep 'Eggers2'],fig1);
    saveas(fig1,[dirData filesep 'Eggers2.fig']);

    
    %% Right hand side of ODE
    function y = ODE(hP)
        
        h = IntM*hP;        
        y = delta./(h.^2+lambda*h) + DDy*hP - hP + theta;
        
        % Boundary conditions 
        y(1)   = hP(1)-1;
        y(end) = hP(end)-theta;
        
    end

    %% f1(x)
    function [y,coeff0_0,coeff1_0] = f1(x,th,thAP)
       if(nargin < 3)
           th   = theta;
           thAP = thetaAP;
       end
       hh0   = th*x+(th-thAP)*(exp(-x)-1);
       y     = -1./(hh0.^2);
       
       % second coefficient for expansion around x = 0
       coeff0_0 = -1/thAP^2;
       coeff1_0 = (th-thAP)/thAP^3;
    end

    function [y,coeff0_0] = f2(x)
        hh0   = theta*x+(theta-thetaAP)*(exp(-x)-1);
        IP    = HIS.InterpolationMatrix_Pointwise(x);        
        y     = 2*(IP*h1)./(hh0.^3);
        
        coeff0_0 = 2*(Cout-1)/thetaAP^5;
    end

    %% Auxiliary functions
    function thAP = ComputeThetaAP()
        thAP = fsolve(@eqn1,theta);
    end

    function y = eqn1(thAP)
        y = 1 - 3*delta*(log(lambda)+ComputeCout(theta,thAP)- Cin )-thAP^3;
    end

    function C = ComputeCout(th,th_ap)
        
        %analytical expression for th=th_ap = thetaAP
%        C         = exp(-1)-expint(1) -1;         
%        for nn = 2:100
%            C = C - (-1)^(nn-1)/(factorial(nn)*(nn-1));
%        end               
        
        ytilde     = 1 + y;
        integrand1 = (1-exp(-ytilde)).*f1(ytilde,th,th_ap);         
        
        ytilde                 = SL.Pts.y;
        [fh,coeff0_0,coeff1_0] = f1(ytilde,th,th_ap);
        integrand2             = (1-exp(-ytilde)).*fh+ 1./((th_ap^2)*ytilde);
        integrand2(1)          = coeff1_0-coeff0_0/2;
                
        C   = (th_ap^2)*(Int*integrand1 + IntL*integrand2) + 1;

    end

    function st = ComputeNextOrder(hhh,st)
        
        f = (st.f);
   
        %% Outer region - 1st term
        % 
        % $$-\frac{1}{2}\int_x^\infty e^{x-t}f(t) dt$$
        % 
        T1 = zeros(N,1);
        for j = 1:N
            t     = y;
            T1(j) = -Int*(exp(-t).*f(y(j)+t))/2;
        end

        %% Outer region - 2nd term
        % 
        % $$-\frac{1}{2}\int_1^x e^{t-x} f(t)dt$$
        %        
        T2 = zeros(N,1);
        T4 = zeros(N,1);
        T5 = zeros(N,1);
        for j = 1:(N-1)
            t     = 1 + SL.Pts.y*(y(j)-1);
            t(1)  = 1;
            ft    = f(t);
            T2(j) = (1-y(j))/2*IntL*(exp(t-y(j)).*ft);
            T4(j) = (1-y(j))*IntL*ft;
            T5(j) = (y(j)-1)*IntL*((cosh(t-y(j))-1).*ft);
        end        
        T2(N) = 0;
        T4(N) = -Int*f(1+y);        
    
        %% Outer region - 3rd term: -K*exp(-x)   
        % 
        % $$K = \int_1^\infty (\frac{e^{-t}}{2} -1) f(t) dt + \int_0^1 (\cosh(t)-1) f(t)dt$$
        %        
        t  = y + 1;    
        ft = f(t);
        T6 = -exp(y)/2*(Int*(exp(-t).*ft));
        K2 = Int*ft;
                
        K  = Int*((exp(-t)/2-1).*ft);
        
        tt           = SL.Pts.y;
        [f_h,c0]     = f(tt);
        
        if(strcmp(func2str(f),'LubricationForcedWetting/f1'))                    
            integrand    = (cosh(tt)-1).*f_h;
            integrand(1) = c0/2;            
        elseif(strcmp(func2str(f),'LubricationForcedWetting/f2'))
            integrand    = (cosh(tt)-1).*f_h+2*log(tt)/thetaAP^5;
            integrand(1) = c0/2;  
            
            K   = K + 2/thetaAP^5;
        end
        K            = K + IntL*integrand;
            
        m1       = (y<1);
        hN(~m1)  = T1(~m1) - T2(~m1) + T4(~m1) + K*exp(-y(~m1)) + K2;
        hN(m1)   = T6(m1) + T5(m1) + K*exp(-y(m1)) + K2;
        hN(end)  = 0;
        hN(1)    = 0;
        hN       = hN';
        
        hNP      = T1 + T2 - K*exp(-y);        
        
        st = v2struct(hN,hNP);
    end

end

