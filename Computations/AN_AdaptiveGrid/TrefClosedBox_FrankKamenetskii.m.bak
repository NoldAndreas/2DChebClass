function TrefClosedBox_FrankKamenetskii()
%************************************************
%
% Solves 
%   (DYN 1) du/dt = Lap(u)+e^u
%   (BC)    u     = 0
%************************************************

    disp(' ** Adaptive Grid Frank Kamenetskii Equation **');
    close all;
    
    %************************************************
    %****************  Preprocess  ****************
    %************************************************        
    N1          = 3; 
    N2          = 56;
    PlotArea    = struct('y1Min',-1,'y1Max',1,'N1',50,...
                         'y2Min',-1,'y2Max',1,'N2',50);
                     
	y2_1D    = (-1:0.01:1)';     
              
    TB                        = BoxTrefSpectralSpectral(N1,N2);
    [Pts,Diff,Int,Ind,Interp] = TB.ComputeAll(PlotArea);        
    Interp1D = TB.InterpolationMatrix_Pointwise(ones(size(y2_1D)),y2_1D);
    Interp0  = TB.InterpolationMatrix_Pointwise(1,0);
    
    %***********************************************
    %*************** Initial Condition *************
    %***********************************************
    T_Adapt    = 0;
    t_record   = [3.544,3.54466,3.5446645,3.544664598];
    params.Dt  = 0.2;
    params.t_n = 0;
    u_n = zeros(size(Pts.y1_kv));
      
    %TB.doPlots(u_n,'SC');

    %****************************************************************
    %****************  Compute time-dependent solution   ************
    %****************************************************************
	mM               = ones(N1*N2,1);
    mM(Ind.bound)    = 0;
    opts             = odeset('RelTol',10^-13,'AbsTol',10^-13,'Mass',diag(mM));           
    
    figure('Color','white','Position',[0 0 750 750]);
    k = 1;
    xNear0 = (-1:0.002:1)';
    
    for i = 1:4000        
        
        st = ComputeStep(params,u_n);
        
        u_n        = st.u_n;
        params.t_n = st.t_n;
        
        if(st.t_n > T_Adapt)
            TB         = st.TB;        
            Diff       = TB.Diff;
            Interp1D   = st.Interp1D;
            Interp0    = st.Interp0;
            params.Dt  = 0.02/max(abs(ft(params.t_n,u_n)));
            params.Dt  = min(params.Dt,t_record(1));
        end        
        
        %[u_n,t_n] = odeSolver(u_n,t_n,Dt);
        %[u_n,t_n] = EulerForward(u_n,t_n,Dt);
        hold off;        
        
        %***************************************************
        %Do Plots
        subplot(2,2,1); hold off;
        mark = (TB.Pts.y1_kv == 1);        
        plot(y2_1D,Interp1D*u_n,'k','linewidth',2); hold on;
        h = plot(TB.Pts.y2_kv(mark),u_n(mark),'o');
        if(st.t_n > 3.3)
            plot(xNear0,AnalyticalSolution(xNear0,st.t_n),'b--');
        end
        set(h,'MarkerEdgeColor','k','MarkerFaceColor','g');
        set(gca,'fontsize',20);                        
        set(gca,'linewidth',1.5);   
        ylim([0 20]);        
        
        title(['t = ',num2str(params.t_n,10)]);
        
        subplot(2,2,2); hold on;
        u0 = Interp0*u_n;
        plot(params.t_n,u0,'o');
        ylabel('u(y_2=0)');
        xlabel('t');
        disp(['t = ',num2str(params.t_n),' , u(0) = ',num2str(u0)]);
        
        %Adapt Grid
        if(params.t_n>=T_Adapt)
            subplot(2,2,3);
            hold on;
            plot(params.t_n,TB.ep_x1(1),'o');                        
                        
            disp(['Dt = ',num2str(params.Dt)]);          
        end
        pause(0.05);
        %***************************************************
        
        if(params.t_n >= t_record(1))
            break;
        end
    end          
    
    function st = ComputeStep(params,u_n)        
        [st.u_n,st.t_n] = odeSolver(u_n,params.t_n,params.Dt);
        if(st.t_n > T_Adapt)
            st.u_n      = TB.UpdatePadeValues(st.u_n);                       
            st.TB       = TB;
            st.Interp1D = TB.InterpolationMatrix_Pointwise(ones(size(y2_1D)),y2_1D);        
            st.Interp0  = TB.InterpolationMatrix_Pointwise(1,0);
        end                
    end
    
    function dudt = ft(t,u)        
        dudt            = Diff.DDy2*u+exp(u);        
        
        dudt(Ind.bound) = u(Ind.bound);        
        dudt(Ind.right) = Ind.normalRight*Diff.grad*u;
        dudt(Ind.left)  = Ind.normalLeft*Diff.grad*u;        
    end

    function u = AnalyticalSolution(x,t)
        T   = 3.544664598;
        xi  = - log(T-t);
        eta = x/sqrt(4*xi*(T-t));
        u   = xi - log(1+eta.^2) - 5/2*(log(xi)/xi).*(eta.^2)./(1+eta.^2);
    end

    function [u_nP1,t_nP1]=odeSolver(u_nM1,t_nM1,Dt)         
%         times = t_nm1;
%         for ii = 1:length(t_record)
%             if((t_record(ii) < t_nM1+Dt) && (t_record(ii) > t_nM1))
%                 times(ii+1) = t_record(ii);
%             end
%         end
        times            = [t_nM1,t_nM1+Dt];
        
        [outTimes,u_t]   = ode15s(@ft,times,u_nM1,opts);
        t_nP1            = outTimes(end);
        u_nP1            = u_t(end,:)';
    end

end