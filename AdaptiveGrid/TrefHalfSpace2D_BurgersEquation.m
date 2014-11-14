function TrefHalfSpace2D_BurgersEquation()
%************************************************
%
% Solves 
%   (DYN 1) du/dt = - u du/dy + nu*Lap(u)
%   (BC)    u     = 0
%************************************************

    disp(' ** Adaptive Grid Burgers Equation Half Space **');
    close all;    
    computeAll = false;
    
    %************************************************
    %****************  Preprocess  ****************
    %*********************************************    
    nu          = 0.001;
    params.nu   = nu;    
       
    N1          = 20;
    N2          = 50;
    params.PhysSpace   = struct('N',[N1,N2],'L1',2,'L2',5);
    PlotArea    = struct('y1Min',-4,'y1Max',4,'N1',50,...
                         'y2Min',0,'y2Max',10,'N2',50);
              
    THS                       = TrefHalfSpace(params.PhysSpace);
    [Pts,Diff,Int,Ind,Interp] = THS.ComputeAll(PlotArea);    
    
    %***********************************************
    %*************** Initial Condition *************
    %***********************************************
    
    params.Dt    = 0.1;
    params.t_n   = 0;    
    params.Identifier = 't_n';
    
    
    %A            = 1+0.5*exp(-(Pts.y1_kv/1).^2);
    A            = 1+0.2*tanh(Pts.y1_kv);
    xh           = (1+Pts.x2_kv)/2;
    u_n          = A.*(sin(2*pi*xh) + 0.5*sin(pi*xh));              
      
    
    THS.PlotXGrid();
    
    THS.plot(u_n,'SC');
    
    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
	mM             = ones(N1*N2,1);
    mM(Ind.bottom) = 0;
    opts           = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));       
                    
    u_n = DataStorage('TrefHalfSpace2D_BurgersEquation',@ComputeStep,params,u_n,computeAll);        
    params.t_n = params.t_n + params.Dt;
    
     for i = 1:200
                
        %Adapt Grid                    
        [u_n,Pts,Diff,Int,Ind,Interp] = THS.UpdatePadeValues(u_n,PlotArea);
        
        hold off;        
        THS.plot(u_n,'SC');
        THS.PlotLineOfPoles(u_n,PlotArea);
        view([1 1 1]);

        %THS.plot(u_n,'contour');
        %THS.PlotLineOfPoles();
        title(['t = ',num2str(params.t_n)]);
        pause(0.05);                
        
        params.Dt  = roundsd(0.3/max(abs(ft(params.t_n,u_n))),3);
        disp(['Dt = ',num2str(params.Dt)]);                              
                
        u_n        = DataStorage('TrefHalfSpace2D_BurgersEquation',@ComputeStep,params,u_n,computeAll);    
        params.t_n = params.t_n + params.Dt;
        
    end                  
    
	function dudt = ft(t,u)                       
        dudt             = nu*Diff.Lap*u-(diag(u)*Diff.Dy2)*u;
        dudt(Ind.bottom) = u(Ind.bottom);
    end
    
    function [u_nP1,t_nP1] = odeSolver(u_nM1,t_nM1,Dt)                
        [outTimes,u_t]     = ode15s(@ft,[t_nM1,t_nM1+Dt],u_nM1,opts);
        t_nP1              = outTimes(end);
        u_nP1              = u_t(end,:)';
    end

    %***********************************
    %******* Auxiliary Functions *******
    %***********************************
    function u_n = ComputeStep(params,u_n)        
        [u_n] = odeSolver(u_n,params.t_n,params.Dt);
    end

end