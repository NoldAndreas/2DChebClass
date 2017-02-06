function DiffusionPeriodicSlit()

    disp('** DiffusionPeriodicSlit **');
    
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
    N1          = 20;       N2      = 30;
    shape.N = [N1;N2];  
              
    abox               = PeriodicBox(shape);
    [Pts,Diff,Int,Ind] = abox.ComputeAll();    
    Interp             = abox.ComputeInterpolationMatrix(...
                                     (-1:0.02:1)',(0:0.02:1)',true,true);                                            
                                               
    rho_ic = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
         
    n1     = length(Pts.x1);
    n2     = length(Pts.x2);

    abox.plot(rho_ic);    

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************        
    mM  = ones(n1*n2,1);
    mM(Ind.left) = 0;
    mM(Ind.right) = 0;
    %opts = odeset('OutputFcn',@PlotFcn,'RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    opts = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.1:2],rho_ic,opts);
    
    for i=1:length(outTimes)
        rho = Rho_t(i,:)';
           
        z = Interp.InterPol*rho;
        
        subplot(2,1,1)
        abox.plot(rho);
        title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);       
        zlim([0 max(rho_ic)]);
        
        subplot(2,1,2)
        h = max( rho -  AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,outTimes(i)));
        hI = max( z -  AnalyticalSolution(Interp.pts1,Interp.pts2,outTimes(i)));
        hold on;
        plot(outTimes(i),h,'o');
        hold on;        
        plot(outTimes(i),hI,'r');
        
        pause(0.02);
        
    end
        
    
    function dydt = Lap(t,rho)
        dydt            = Diff.Lap*rho;
        dydt(Ind.left) = Ind.normalLeft*(Diff.grad*rho);
        dydt(Ind.right) = Ind.normalRight*(Diff.grad*rho);
    end
                   
    function y = AnalyticalSolution(y1,y2,t)
        y = cos(pi*y1/L).*cos(2*pi*y2/L)*exp(-5*t*(pi/L)^2);
    end

    function stat = PlotFcn(t,rho,flag)
        if(~ isscalar(t))
            t = 0;
        end        
        
        stat = 0;
 
        if isempty(flag)
            
            rho_full = rho;
            
            doPlots_IP(Interp,rho_full);      
            zlim([0 max(rho_ic)]);
            drawnow;
        end
    end


end