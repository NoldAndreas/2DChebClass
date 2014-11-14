function DiffusionWedge()

    disp('** Diffusion Wedge **');

	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    L  = 3;
    N1 = 40; 
    N2 = 20;
    half_wedge_angle    = pi/4;
    Geometry = struct('R_in',0,'LR',L,...
                      'th1',-half_wedge_angle,...
                      'th2',half_wedge_angle,'N',[N1;N2]);
    
    WDG                = InfWedge(Geometry);
    [Pts,Diff,Int,Ind] = WDG.ComputeAll();
    Interp = WDG.ComputeInterpolationMatrix((-1:0.02:0.8)',(-1:0.02:1)',true,true);
              
    [rho_ic,Lap_ic]       = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
    [rho_ic_IP,Lap_ic_IP] = AnalyticalSolution(Interp.pts1,Interp.pts2,0);    
    
    figure
    WDG.plot(rho_ic);
    figure
    subplot(2,1,1); WDG.plot(Lap_ic);
    subplot(2,1,2); WDG.plot(Interp.InterPol*Lap_ic-Lap_ic_IP);    
    
    errIP = Interp.InterPol*rho_ic - rho_ic_IP;
    display(['error Interpolation: ' , num2str(max(abs(errIP)))]);
                   
    errLap = Diff.Lap(~Ind.bound,:)*rho_ic - Lap_ic(~Ind.bound);
    display(['error laplace: ' , num2str(max(abs(errLap)))]);
        
    %Compute matrix of Integration over subspace            
    subGeometry   = struct('R_in',0,'R_out',L,...
                      'th1',-half_wedge_angle,...
                      'th2',half_wedge_angle,'N',[20;20]);    
    subShape      = Wedge(subGeometry);
    IP            = WDG.SubShapePts(subShape.Pts);
    Int_of_path   = subShape.IntFluxThroughDomain(100)*blkdiag(IP,IP);
    Int_SubOnFull = subShape.ComputeIntegrationVector()*IP;
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;    
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.1:5],[0;rho_ic],opts);
    
    %************************************************
    %****************  Postprocess  ****************
    %************************************************
    
    figure
    for i=1:length(outTimes)
        rho     = Rho_t(i,2:end)';        
        accFlux = Rho_t(i,1);
        t       = outTimes(i);
        
        z     = Interp.InterPol*rho;
        z_ana = AnalyticalSolution(Interp.pts1,Interp.pts2,outTimes(i));
        
        subplot(2,2,1)
        WDG.plot(rho);        
        title(['Interpolation of Solution at t = ', num2str(t)]);       
        zlim([0 max(rho_ic)]);
        
        subplot(2,2,2)
        hold on;
        h = max( abs(rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,outTimes(i))));
        hI = max( abs(z -  z_ana));
        hold on;
        plot(outTimes(i),h,'o',outTimes(i),hI,'r');
        plot(outTimes(i),accFlux + Int_SubOnFull*(rho-rho_ic),'ob'); hold on;
        plot(outTimes(i),Int_SubOnFull*(rho-rho_ic),'or');  hold on;
        plot(outTimes(i),accFlux,'og');
        title(['Mass in the system. Deviation of analytical Sol: ' , num2str(hI)]);
        
        subplot(2,2,3)
        flux = -Diff.grad*rho;
        WDG.plotFlux(flux); hold on;        
        subShape.PlotBorders();
        
        subplot(2,2,4);   
        hold off;
        subShape.borderRight.PlotValueOnPath(blkdiag(IP,IP)*flux,'normal');
        [h1s,h2s,h3s,h4s,fluxInt] = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,t,L);
        title(['Error of Flux (compare to Ana. Sol.): ' , ...
                 num2str(Int_of_path*flux - fluxInt)]);        
         xlabel('\theta in deg');
        
        pause(0.02);
        
    end
        
    
    function dydt = Lap(t,rho)        
        rho              = rho(2:end);    
        
        dydt             = Diff.Lap*rho;
        
        flux             = -Diff.grad*rho;
        
        %no flux at origin
        dydt(Ind.left)  = Ind.normalLeft*flux;
        %zero value at infinity
        dydt(Ind.right) = rho(Ind.right);
        %no flux at top and bottom wedge
        dydt(Ind.top & ~Ind.left & ~Ind.right)   = Ind.normalTop(2:end-1,:)*flux;
        dydt(Ind.bottom & ~Ind.left & ~Ind.right)   = Ind.normalBottom(2:end-1,:)*flux;                       
        
        dydt = [Int_of_path*flux;dydt];
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