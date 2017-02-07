function DiffusionWedge()

    disp('** Diffusion in an infinite wedge **');
    AddPaths();
    close all;
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************    
    half_wedge_angle    = pi/4;
    Geometry = struct('R_in',0,'LR',3,...
                      'th1',-half_wedge_angle,...
                      'th2',half_wedge_angle,'N',[40;20]);
    
    WDG              = InfWedge(Geometry);
    [Pts,Diff,~,Ind] = WDG.ComputeAll();
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM               = ones(WDG.M,1);
    mM(Ind.bound)    = 0;    
    opts             = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [rho_ic,Lap_ic]  = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.1:5],rho_ic,opts);
    
    %******************************************
    %Compare analytical vs. numerical solution
    Interp          = WDG.ComputeInterpolationMatrix((-1:0.02:0.25)',(-1:0.02:1)',true,true);        
    [rho_ic_IP,~]   = AnalyticalSolution(Interp.pts1,Interp.pts2,0);    
           
    errIP = Interp.InterPol*rho_ic - rho_ic_IP;
    display(['Error of interpolation: ' , num2str(max(abs(errIP)))]);
                   
    errLap = Diff.Lap(~Ind.bound,:)*rho_ic - Lap_ic(~Ind.bound);
    display(['Error of Laplace operator: ' , num2str(max(abs(errLap)))]);
    
    %****************************************
    %****************  Plot  ****************
    %****************************************
    
    fig1 = figure;
    for i=1:length(outTimes)
        rho     = Rho_t(i,:)';        
        t       = outTimes(i);
        
        %*** Compute error wrt analytical solution ***
        z             = Interp.InterPol*rho;
        z_ana         = AnalyticalSolution(Interp.pts1,Interp.pts2,outTimes(i));        
        error_anal(i) = max( abs(z -  z_ana));
        
        set(fig1,'Name',['Interpolation of Solution at t = ', num2str(t)]);
        set(fig1,'NumberTitle','off');
        %*** Solution at time t ***
        subplot(1,2,1)
        WDG.plot(rho);                
        zlim([0 max(rho_ic)]);        
        
        %*** Flux at time t ***
        subplot(1,2,2); hold off;
        WDG.plotFlux(-WDG.Diff.grad*rho); 
        
        pause(0.05);        
    end
    
    figure;
    plot(outTimes,error_anal);    
    xlabel('t');
    ylabel('Error');
        
    
    function dydt = Lap(t,rho)        
        
        dydt             = Diff.Lap*rho;        
        flux             = -Diff.grad*rho;
        
        %no flux at origin
        dydt(Ind.left)  = Ind.normalLeft*flux;
        %zero value at infinity
        dydt(Ind.right) = rho(Ind.right);
        %no flux at top and bottom wedge
        dydt(Ind.top & ~Ind.left & ~Ind.right)    = Ind.normalTop(2:end-1,:)*flux;
        dydt(Ind.bottom & ~Ind.left & ~Ind.right) = Ind.normalBottom(2:end-1,:)*flux;                                       
    end
                   
    function [y,Lap,dydr,dydtheta,fluxInt] = AnalyticalSolution(r,theta,t,Lf)
        t0 = -1;
        
        x0 = 0;
        
        x = r.*cos(theta);
        y = r.*sin(theta);
        
        r  = sqrt((x-x0).^2 + y.^2);
        y  = exp(-r.^2./(4*(t-t0)))./(t-t0);        
        
        Lap = ((r.^2 - 4*(t -t0))/(4*(t-t0)^2)).*y;        
        Lap(r == inf) = 0;
        
        dydr      = -r.*y/(2*(t-t0));
        dydtheta  = zeros(size(r));
        
        if(nargin == 4)
            fluxs   = Lf/(2*(t-t0))*exp(-Lf.^2./(4*(t-t0)))./(t-t0);    
            fluxInt = fluxs*2*half_wedge_angle*Lf;
        end
    end

end