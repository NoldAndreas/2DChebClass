function DiffusionClosedBoxTest()
%************************************************
%DiffusionClosedBox()
%
% Solves 
%   (DYN 1) dx/dt     = Lap(x)
%   (BC 1)  n*grad(x) = 0 
%************************************************

    disp(' ** DiffusionClosedBox **');

	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************
    L = 2;
    shape.y1Min = 0;	shape.y1Max = L;	
    shape.y2Min = 0;	shape.y2Max = L;	    
    %N1          = 21;       N2      = 20;
    N1          = 50;       N2      = 50;
    shape.N     = [N1;N2];  
              
    abox               = Box(shape);
    [Pts,Diff,Int,Ind] = abox.ComputeAll();    
    Interp             = abox.ComputeInterpolationMatrix(...
                                    (-1:0.02:1)',(-1:0.02:1)',true,true);                      

    %rho_ic             = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);               
    rho_ic             = InitialCondition(Pts.y1_kv,Pts.y2_kv);               
    
    %Compute matrix of Integration over subspace
    shapeSub      = struct('y1Min',0,'y1Max',L/2,'y2Min',0,'y2Max',L,'N',[20,20]);
    subBox        = Box(shapeSub);
    IP            = abox.SubShapePts(subBox.Pts);
    Int_SubOnFull = subBox.ComputeIntegrationVector()*IP;
        
    Int_of_path   = subBox.IntFluxThroughDomain(100)*blkdiag(IP,IP);
    
    h              = subBox.borderRight.InterpOntoBorder*IP;
    RightBorderPts = subBox.borderRight.normal*blkdiag(h,h);
    
    interpBR = subBox.borderRight.ComputeInterpolationMatrix((-1:0.02:1)'); %BR = Border Right
    IntBR    = subBox.borderRight.IntNormal*blkdiag(IP,IP);

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM               = ones(N1*N2,1);
    mM(Ind.bound)    = 0;
    opts             = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag([1;mM]));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.05:1],[0;rho_ic],opts);

    %************************************************
    %****************  Postprocess  ****************
    %************************************************        
    
    for i=1:length(outTimes)
        rho     = Rho_t(i,2:end)';  
        accFlux = Rho_t(i,1);
        t       = outTimes(i);
        
        subplot(2,2,1)
        abox.plot(rho);
        title(['Interpolation of Solution at t = ', num2str(t)]);       
        zlim([0 max(rho_ic)]);
                
        
        %Check mass conservation
        subplot(2,2,2)   
        h  = abs(max( rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,t)));
        hI = abs(max( Interp.InterPol*rho - ...
            AnalyticalSolution(Interp.pts1,Interp.pts2,t)));        
        hold on;
        plot(t,accFlux + Int_SubOnFull*(rho-rho_ic),'ob'); hold on;
        plot(t,Int_SubOnFull*(rho-rho_ic),'xr');  hold on;
        plot(t,accFlux,'+g'); hold on;
        plot(t,Int*(rho - rho_ic),'om');
        xlabel('t');
        legend('Error in mass of subsystem',...
                'Mass gain in subsystem',...
                '(Out) flux from subsystem',...
                'Mass gain in full system');
        title(['Mass in the full and subsystem. Error to Anal. Sol.: '...
            num2str(max(h,hI))]);                
        
        subplot(2,2,3);
        fl   = GetFlux(rho,outTimes(i));
        abox.plotFlux(fl); hold on;%NormQuiverPlot(Pts,fl);
        subBox.PlotBorders();
        
        subplot(2,2,4);   
        hold off;
        subBox.borderRight.PlotValueOnPath(blkdiag(IP,IP)*fl,'normal');
        [z,dz_dy1,dz_dy2] = AnalyticalSolution(subBox.borderRight.Pts.y1_kv,...
                                           subBox.borderRight.Pts.y2_kv,t);
        plot(subBox.borderRight.Pts.t,-dz_dy1,'g');
        title(['Flux through the path is: ' , num2str(IntBR*fl)]);
        xlabel('t');
        
        pause(0.01);
        pause
        drawnow;
    end
    
    
    function dydt = Lap(t,rho)
        rho             = rho(2:end);
        dydt            = Diff.Lap*rho;
        flux            = GetFlux(rho,t);
        dydt(Ind.bound) = Ind.normal*flux;
        dydt            = [Int_of_path*flux;dydt];
    end

    function flux = GetFlux(rho,t)        
        flux  = -Diff.grad*rho;                        
    end

                   
    function [z,dz_dy1,dz_dy2] = AnalyticalSolution(y1,y2,t)       
        z      = cos(pi*y1/L).*cos(pi*y2/L)*exp(-2*t*(pi/L)^2);
        dz_dy1 = -pi/L*sin(pi*y1/L).*cos(pi*y2/L)*exp(-2*t*(pi/L)^2);
        dz_dy2 = -pi/L*cos(pi*y1/L).*sin(pi*y2/L)*exp(-2*t*(pi/L)^2);
        %z  = 1-(y1/L).^2 - (y2/L).^2;
        %dz_dy1 =  -2*y1/L;
        %dz_dy2 =  -2*y2/L;
    end

    function z = InitialCondition(y1,y2)
        z = zeros(size(y1));
        z(y1>0.5) = 1;
    end

end