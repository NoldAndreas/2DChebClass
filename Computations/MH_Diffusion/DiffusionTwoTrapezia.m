function DiffusionTwoTrapezia()

    disp('** DiffusionTwoTrapezia **');
    
    %---------------------------------------------------------------------%
    % Set up physical space                                               %
    %---------------------------------------------------------------------%

    N1 = 12;
    N2 = 12;
    
    geom.L1 = 1;
    geom.L2 = 3;
    geom.h = 2;
    geom.N = [N1;N2];
    
    trap  = Trapezium2(geom);
    
    geomT.L1 = 2;
    geomT.L2 = 3;
    geomT.h = 1;
    geomT.N = [N1;N2];
    geomT.shift = [0;0.5]; % centre of trapezium
    
    NT = geomT.N(1)*geomT.N(2);
    
    geomB.L1 = 1;
    geomB.L2 = 2;
    geomB.h = 1;
    geomB.N = [N1;N2];
    geomB.shift = [0;-0.5]; % centre of trapezium
    
    NB = geomB.N(1)*geomB.N(2);
    
    trapT = Trapezium2(geomT);
    trapB = Trapezium2(geomB);

    interpGrid = (-1:0.05:1)';
    
    [Pts,Diff,Int,Ind] = trap.ComputeAll();    
    Interp             = trap.ComputeInterpolationMatrix(interpGrid,interpGrid,true,true);
   
    [PtsT,DiffT,IntT,IndT] = trapT.ComputeAll();    
    InterpT             = trapT.ComputeInterpolationMatrix(interpGrid,interpGrid,true,true);
    
    [PtsB,DiffB,IntB,IndB] = trapB.ComputeAll();    
    InterpB             = trapB.ComputeInterpolationMatrix(interpGrid,interpGrid,true,true);
    
    %---------------------------------------------------------------------%
    % Make intitial condition                                             %
    %---------------------------------------------------------------------%
    
    rho_ic = InitialCondition(Pts.y1_kv,Pts.y2_kv,0);

    rho_icT = InitialCondition(PtsT.y1_kv,PtsT.y2_kv,0);
    
    rho_icB = InitialCondition(PtsB.y1_kv,PtsB.y2_kv,0);
    
    rho_icTB = [rho_icT; rho_icB];

    m_ic = Int*rho_ic;
    
    m_icT = IntT*rho_icT;
    m_icB = IntB*rho_icB;
    
    m_iCTB = m_icT + m_icB;
         
%     figure
%     subplot(1,2,1)
%     trap.plot(rho_ic);
%     axis tight
%     
%     subplot(1,2,2)
%     trapT.plot(rho_icT);
%     hold on
%     trapB.plot(rho_icB);
%     axis tight
%     return
    
    
    %---------------------------------------------------------------------%
    % Solve the ODE                                                       %
    %---------------------------------------------------------------------%
        
    mM              = ones(size(Pts.y1_kv));
    mM(Ind.bound) = 0;    
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.001:0.1],rho_ic,opts);

    
    
    mMT              = ones(size(PtsT.y1_kv));
    mMT(IndT.bound) = 0;
    mMB              = ones(size(PtsB.y1_kv));
    mMB(IndB.bound) = 0;

    mMTB = [mMT;mMB];
    optsTB = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mMTB));
    [outTimesTB,Rho_tTB] = ode15s(@LapTB,[0:0.001:0.1],rho_icTB,optsTB);

        
    %---------------------------------------------------------------------%
    % Plot the solutions                                                  %
    %---------------------------------------------------------------------%
   
     figure
    
    for i = 1:length(outTimes)
        
        
        rhot = Rho_t(i,:)';
        
        rhotTB = Rho_tTB(i,:)';
        
        rhotT = rhotTB(1:NT);
        rhotB = rhotTB(NT+1:NT+NB);
        
       
        subplot(1,2,1);
        hold off
        trapT.plot(rhotT);
        hold on
        trapB.plot(rhotB);
        hold on
        trap.plot(rhot);
        axis tight
        
        
        subplot(1,2,2);

        hold off
        plot(Pts.y1_kv(Ind.bottom),rhot(Ind.bottom),'b-o')
        hold on
        plot(PtsB.y1_kv(IndB.bottom),rhotB(IndB.bottom),'r--x')


        plot(Pts.y1_kv(Ind.top),rhot(Ind.top),'m-o')
        hold on
        plot(PtsT.y1_kv(IndT.top),rhotT(IndT.top),'g--x')

        
        
%         subplot(1,2,2)
%         
%         trap.plot(rhot);
%         axis tight
             
%         xlim([0;2]);
%         ylim([0;1]);
        
       pause(0.05);
        
       %pause
    end
    
    %---------------------------------------------------------------------%
    % Right hand side of ODE                                              %
    %---------------------------------------------------------------------%

    function dydt = Lap(t,rho)
        dydt            = Diff.Lap*rho;
        dydt(Ind.bound)  = Ind.normal*(Diff.grad*rho); % compute flux directly
    end

    function dydt = LapTB(t,rho)
        
        Tmask = 1:NT;
        Bmask = NT + 1 : NT+NB;
        
        rhoT    = rho(Tmask);
        rhoB    = rho(Bmask);
        
        dydtT = DiffT.Lap*rhoT;
        dydtB = DiffB.Lap*rhoB;

        fluxTTop = IndT.normalTop*(DiffT.grad*rhoT);
        fluxTLeft = IndT.normalLeft*(DiffT.grad*rhoT);
        fluxTRight = IndT.normalRight*(DiffT.grad*rhoT);
        
        fluxBBottom = IndB.normalBottom*(DiffB.grad*rhoB);
        fluxBLeft = IndB.normalLeft*(DiffB.grad*rhoB);
        fluxBRight = IndB.normalRight*(DiffB.grad*rhoB);

        dydtT(IndT.top) = fluxTTop;
        dydtT(IndT.left) = fluxTLeft;
        dydtT(IndT.right) = fluxTRight;
        
        dydtB(IndB.bottom) = fluxBBottom;
        dydtB(IndB.left) = fluxBLeft;
        dydtB(IndB.right) = fluxBRight;
        
        dydtT(IndT.bottom) = rhoT(IndT.bottom) - rhoB(IndB.top);
        dydtB(IndB.top) = IndT.normalBottom*(DiffT.grad*rhoT) + IndB.normalTop*(DiffB.grad*rhoB);
               
        dydt = [dydtT; dydtB];
                
    end


    function y = InitialCondition(y1,y2,t)
        
        %y = sin(y1)+1;
        
        y = exp(-(y1+0.1).^2/0.5-(y2-0.25).^2)/0.5;
        %y = cos(y1*2*pi).*cos(y2*2*pi) + 1;
        
        %yb = exp(-(y1-1).^2-(y2-0.5).^2);
        %yt = exp(-(y1-1).^2-(y2-0.5).^2);
        %y = [yb; yt];
    end

end