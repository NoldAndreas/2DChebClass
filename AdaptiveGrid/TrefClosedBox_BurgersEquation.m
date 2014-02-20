function TrefClosedBox_BurgersEquation()
%************************************************
%
% Solves 
%   (DYN 1) du/dt = Lap(u)+e^u
%   (BC)    u     = 0
%************************************************

    disp(' ** Adaptive Grid Burgers Equation **');
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************    
    %K           = 1;
    nu          = 0.001;
    
    N1          = 3; 
    N2          = 50;
    PlotArea    = struct('y1Min',-1,'y1Max',1,'N1',50,...
                         'y2Min',-1,'y2Max',1,'N2',50);
              
    TB                        = BoxTrefSpectralSpectral(N1,N2);
    [Pts,Diff,Int,Ind,Interp] = TB.ComputeAll(PlotArea);
    y2_1D      = (-1:0.01:1)';     
    Interp1D   = TB.InterpolationMatrix_Pointwise(ones(size(y2_1D)),y2_1D);
    Interp1D_C = TB.ComputeInterpolationMatrix(1,y2_1D,true);
    
    mark1 = (1:N1*N2);
    mark2 = (1+N1*N2:2*N1*N2);
    
    %***********************************************
    %*************** Initial Condition *************
    %***********************************************
       
    Dt           = 0.05;
    t_n          = 0;    
    xh           = (1+Pts.y2_kv)/2;
    u_n          = (sin(2*pi*xh) + 0.5*sin(pi*xh));    
          %(sin(2*pi*Pts.y1_kv) + 0.5*sin(pi*Pts.y1_kv)).*...
          %(sin(2*pi*Pts.y2_kv) + 0.5*sin(pi*Pts.y2_kv));      	

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
	mM            = ones(N1*N2,1);
    mM(Ind.bound) = 0;
    opts          = odeset('RelTol',10^-8,'AbsTol',10^-8,'Mass',diag(mM));   
    
    figure('Color','white','Position',[0 0 750 750]);
    
    for i = 1:200
        [u_n,t_n] = odeSolver(u_n,t_n,Dt);
        %[u_n,t_n] = EulerForward(u_n,t_n,Dt);
        hold off;
        
        %***************************************************
        %Do Plot
        mark = (TB.Pts.y1_kv == 1);        
        
        subplot(1,2,1); hold off;
        plot(y2_1D,Interp1D*u_n,'k','linewidth',2); hold on;
        h = plot(TB.Pts.y2_kv(mark),u_n(mark),'o');
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
        set(gca,'fontsize',20);                        
        set(gca,'linewidth',1.5);           
        title(['t = ',num2str(t_n)]);
        
        subplot(1,2,2);  hold off;      
        plot(y2_1D,Interp1D_C.InterPol*u_n,'k','linewidth',2); hold on;
        h = plot(TB.Pts.x2_kv(mark),u_n(mark),'o');
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
        set(gca,'fontsize',20);                        
        set(gca,'linewidth',1.5);           
        title(['t = ',num2str(t_n)]);
        
        %***************************************************                
        
        %Adapt Grid
        if(t_n>=0.1)
            %Dt = 0.2;
            [u_n,Pts,Diff,Int,Ind,Interp] = TB.UpdatePadeValues(u_n,PlotArea);
            y2_1D      = (-1:0.01:1)';     
            Interp1D   = TB.InterpolationMatrix_Pointwise(ones(size(y2_1D)),y2_1D);
            Interp1D_C = TB.ComputeInterpolationMatrix(1,y2_1D,true);
            
            %TB.PlotLineOfPoles(u_n(mark2));
            Dt  = 0.3/max(abs(ft(t_n,u_n)));
            disp(['t_n = ',num2str(t_n),'   ,    Dt = ',num2str(Dt)]);          
        end
        pause(0.05);
    end          
    
    function dudt = ft(t,u)       
        dudt   = nu*Diff.Lap*u-diag(u)*Diff.Dy2*u;

        dudt(Ind.bound)  = u(Ind.bound);
        dudt(Ind.right)  = Ind.normalRight*Diff.grad*u;
        dudt(Ind.left)   = Ind.normalLeft*Diff.grad*u;
    end

    function [u_nP1,t_nP1]=odeSolver(u_nM1,t_nM1,Dt)                
        [outTimes,u_t]   = ode15s(@ft,[t_nM1,t_nM1+Dt],u_nM1,opts);
        t_nP1            = outTimes(end);
        u_nP1            = u_t(end,:)';
    end

end