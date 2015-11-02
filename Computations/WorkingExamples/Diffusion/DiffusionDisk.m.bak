function DiffusionDisk()

    disp('** DiffusionPolarDisk **');

	if(length(dbstack) == 1)
        AddPaths();
    end    
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************
    tic
    R  = 3;
    N  = [20,20];
        
    DC                 = Disc(v2struct(R,N));
    [Pts,Diff,Int,Ind] = DC.ComputeAll();    
    Interp  = DC.ComputeInterpolationMatrix((0:0.02:1)',(0:0.02:1)',true,true);
    
    %Check Polar Spectral/Fourier map in 2D:
              
    rho_ic = AnalyticalSolution(Pts.y1_kv,Pts.y2_kv,0);
    rho_IP = AnalyticalSolution(Interp.pts1,Interp.pts2,0);

    subplot(2,1,1); DC.doPlots(rho_ic);
    subplot(2,1,2); DC.doPlots(Interp.InterPol*rho_ic-rho_IP);
    %doPlots_IP_Polar(Interp,rho_ic,rho_IP);

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
        
    mM              = ones(size(Pts.y1_kv));
    mM(Ind.outR) = 0;    
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.2:5],rho_ic,opts);
    toc
    %****************************************************************
    %***********************  Postprocess   *************************
    %****************************************************************
   
    figure
    for i=1:length(outTimes)
           
        DC.doPlots(Rho_t(i,:)');        
        title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);       
        zlim([0 max(rho_ic)]);        
        
        pause(0.02);
        
    end   
    
    function dydt = Lap(t,rho)
        dydt            = Diff.Lap*rho;
        dydt(Ind.outR)  = Ind.normalOutR*(Diff.grad*rho);        
    end
                   
    function y = f(rho)
        y = DdoubleWellFE(rho) - Vext7(Pts.y1_kv,Pts.y2_kv);
        y = y(Ind.outR == 0);
    end

    function dy = DdoubleWellFE(x)
        dy = 2*x;
    end

    function y = AnalyticalSolution(r_o,theta_o,t)
        y = (R^2-r_o.^2).^2;
    end

%     function stat = PlotFcn(t,rho,flag)
%         if(~ isscalar(t))
%             t = 0;
%         end        
%         
%         stat = 0;
%  
%         if isempty(flag)            
%             doPlots_IP_Polar(Interp,rho);      
%             zlim([0 max(rho_ic)]);
%             drawnow;
%         end
%     end   


end