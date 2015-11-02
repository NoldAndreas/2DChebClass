function DiffusionCapillary()

    disp('** Diffusion Capillary **');

	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************
    L = 2;
	%Check Polar Spectral/Fourier map in 2D:              
    shape.y1Min = 0;	shape.y1Max = L;	
    shape.y2Min = 0;	shape.y2Max = L;	    
    N1          = 20;       N2      = 20;
    shape.N = [N1;N2];  
              
    abox               = Box(shape);
    [Pts,Diff,Int,Ind] = abox.ComputeAll();    
    Interp             = abox.ComputeInterpolationMatrix(...
                                    (-1:0.02:1)',(-1:0.02:1)',true,true);                                            
              
    rho_ic = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
         
    n1     = length(Pts.x1);
    n2     = length(Pts.x2);

    abox.doPlots(rho_ic);
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM  = ones(n1*n2,1);
    mM(Ind.bound) = 0;
    %opts = odeset('OutputFcn',@PlotFcn,'RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,0:0.1:3,rho_ic,opts);
    
    %************************************************
    %****************  Postprocess  ****************
    %************************************************    
    for i=1:length(outTimes)
        rho = Rho_t(i,:)';
           
        subplot(2,1,1)
        abox.doPlots(rho);        
        title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);               
        zlim([0 max(max(Rho_t))]);
        
        subplot(2,1,2)
        z = Interp.InterPol*rho;
        h = max( abs(rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,outTimes(i))));
        hI = max( abs(z -  AnalyticalSolution(Interp.pts1,Interp.pts2,outTimes(i))));
        hold on;
        plot(outTimes(i),h,'o');
        hold on;        
        plot(outTimes(i),hI,'r');
        
        pause(0.02);
        
    end
            
    function dydt = Lap(t,rho)
        dydt            = Diff.Lap*rho;
        drho            = Diff.grad*rho;
        dydt(Ind.left)  = Ind.normalLeft*drho;
        dydt(Ind.right) = Ind.normalRight*drho;
        dydt(Ind.bottom)= rho(Ind.bottom);
        dydt(Ind.top)   = Ind.normalTop*drho;
    end
                   
    function y = AnalyticalSolution(y1,y2,t)
        y = cos(pi*y1/L).*sin(pi*y2/(2*L))*exp(-5/4*t*(pi/L)^2);        
    end

end