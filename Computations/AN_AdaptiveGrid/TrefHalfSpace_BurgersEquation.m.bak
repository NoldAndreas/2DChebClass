function TrefHalfSpace_BurgersEquation()
%************************************************
%
% Solves 
%   (DYN 1) du/dt = Lap(u)+e^u
%   (BC)    u     = 0
%************************************************

    disp(' ** Adaptive Grid Burgers Equation Half Space **');
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %*********************************************    
    nu          = 0.001;
       
    N1          = 3; 
    N2          = 40;
    PhysSpace   = struct('N',[N1,N2],'L1',2,'L2',5);
    PlotArea    = struct('y1Min',-4,'y1Max',4,'N1',50,...
                         'y2Min',0,'y2Max',10,'N2',50);
              
    THS                       = TrefHalfSpace(PhysSpace);
    [Pts,Diff,Int,Ind,Interp] = THS.ComputeAll(PlotArea);    
    
    %***********************************************
    %*************** Initial Condition *************
    %***********************************************
       
    Dt           = 0.05;
    t_n          = 0;    
    xh           = (1+Pts.x2_kv)/2;
    u_n          = (sin(2*pi*xh) + 0.5*sin(pi*xh));              
      
    THS.doPlots(u_n,'SC');

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
	mM            = ones(N1*N2,1);
    mM(Ind.bottom) = 0;
    opts          = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));   
    
    for i = 1:200
        [u_n,t_n] = odeSolver(u_n,t_n,Dt);
        %[u_n,t_n] = EulerForward(u_n,t_n,Dt);
        hold off;
        
        THS.doPlots(u_n,'SC');
        title(['t = ',num2str(t_n)]);
        
        %Adapt Grid
        if(t_n>=0.1)
            %Dt = 0.2;
            [u_n,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(u_n,PlotArea);
            THS.PlotLineOfPoles(u_n);
            Dt  = 0.3/max(abs(ft(t_n,u_n)));
            disp(['Dt = ',num2str(Dt)]);              
        end
        view([1 1 1]);
        pause(0.05);
    end          
    
	function dudt = ft(t,u)                       
        dudt             = nu*Diff.Lap*u-(diag(u)*Diff.Dy2)*u;
        dudt(Ind.bottom) = u(Ind.bottom);
    end
    
    function [u_nP1,t_nP1]=odeSolver(u_nM1,t_nM1,Dt)                
        [outTimes,u_t]   = ode15s(@ft,[t_nM1,t_nM1+Dt],u_nM1,opts);
        t_nP1            = outTimes(end);
        u_nP1            = u_t(end,:)';
    end

end