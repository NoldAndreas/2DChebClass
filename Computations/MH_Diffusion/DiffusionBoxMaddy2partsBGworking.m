function [Interp, InterpT, InterpB] = DiffusionBoxMaddy2partsBGworking()

    disp('** DiffusionBox **');
    
    %---------------------------------------------------------------------%
    % Set up physical space                                               %
    %---------------------------------------------------------------------%

    N1 = 12;
    N2 = 12;
    
    geom.y1Min = 0; 
    geom.y1Max = 2 ;   
    geom.y2Min = 0;
    geom.y2Max = 1;
    geom.N = [N1;N2];
    % make the box
    DC                 = Box(geom);
    
    geomT.y1Min = 0;
    geomT.y1Max = 2;
    geomT.y2Min = 0.5;
    geomT.y2Max = 1;
    geomT.N = [N1;N2];
    
    NT = geomT.N(1)*geomT.N(2);
    
    geomB.y1Min = 0;
    geomB.y1Max = 2;
    geomB.y2Min = 0;
    geomB.y2Max = 0.5;
    geomB.N = [N1;N2];
    
    NB = geomB.N(1)*geomB.N(2);
    
    DCT = Box(geomT);
    DCB = Box(geomB);
    
%     figure
%     DC.PlotGrid
%     
%     figure
%     DCT.PlotGrid
%     hold on
%     DCB.PlotGrid
%     return
    
    NInterp1 = 20;
    NInterp2 = 20;

    y1Plot = (linspace(0,2,NInterp1))';
    y2Plot = (linspace(0,1,NInterp2))';
    y2PlotB = y2Plot(1:NInterp2/2);
    y2PlotT = y2Plot(NInterp2/2+1:NInterp2);

    x1Plot = DC.CompSpace1(y1Plot);
    x2Plot = DC.CompSpace2(y2Plot);
 
    x2PlotB = DCB.CompSpace2(y2PlotB);
    x2PlotT = DCT.CompSpace2(y2PlotT);
    
    % compute points and matrices
    [Pts,Diff,Int,Ind] = DC.ComputeAll();    
    Interp  = DC.ComputeInterpolationMatrix(x1Plot,x2Plot,true,true);

    [PtsT,DiffT,IntT,IndT] = DCT.ComputeAll();    
    InterpT  = DCT.ComputeInterpolationMatrix(x1Plot,x2PlotT,true,true);

    [PtsB,DiffB,IntB,IndB] = DCB.ComputeAll();    
    InterpB  = DCB.ComputeInterpolationMatrix(x1Plot,x2PlotB,true,true);
    
    %---------------------------------------------------------------------%
    % Make intitial condition                                             %
    %---------------------------------------------------------------------%

    rho_ic = InitialCondition(Pts.y1_kv,Pts.y2_kv,0);

    rho_icT = InitialCondition(PtsT.y1_kv,PtsT.y2_kv,0);
    
    rho_icB = InitialCondition(PtsB.y1_kv,PtsB.y2_kv,0);
    
    rho_icTB = [rho_icT; rho_icB];

%     figure
%     
% %     subplot(1,2,1)
%     
%     DCT.plot(rho_icT);
%  
%     hold on
% %     subplot(1,2,2)
%    
%     DCB.plot(rho_icB);
% 
%     xlim([0;2]);
%     ylim([0;1]);
%  
%     
%     return
%     
    %Ind.normal*(Diff.grad*rho_ic);
    
%     rho_ic = 1000*rho_ic;


    m_ic = Int*rho_ic;
    
    m_icT = IntT*rho_icT;
    m_icB = IntB*rho_icB;
    
    m_iCTB = m_icT + m_icB;
    %---------------------------------------------------------------------%
    % Solve the ODE                                                       %
    %---------------------------------------------------------------------%
        
    mM              = ones(size(Pts.y1_kv));
    mM(Ind.bound) = 0;    
    opts = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mM));
    [outTimes,Rho_t] = ode15s(@Lap,[0:0.001:0.1],rho_ic,opts);

    
    
    mMT              = ones(size(PtsT.y1_kv));
    mMT(IndT.bound) = 0;
%     mMT(IndT.top) = 0;
%     mMT(IndT.left) = 0;
%     mMT(IndT.right) = 0;
    mMB              = ones(size(PtsB.y1_kv));
    mMB(IndB.bound) = 0;
%     mMB(IndB.bottom) = 0;
%     mMB(IndB.left) = 0;
%     mMB(IndB.right) = 0;

    mMTB = [mMT;mMB];
    optsTB = odeset('RelTol',10^-9,'AbsTol',10^-9,'Mass',diag(mMTB));
    [outTimesTB,Rho_tTB] = ode15s(@LapTB,[0:0.001:0.1],rho_icTB,optsTB);

    
    
    %return
    
    % at each time step
    
    % apply interpolation matrices to the corresponding solutions to
    % interpolate onto the grids we've made
    
    % make tempBT = [reshape(interpolatedRhoB,5,10);reshape(interpolatedRhoT,5,10)];
    % and  temp   = reshape(interpolatedRho,10,10);
   
    % sum(sum(abs(temp-tempBT).^2)) / sum(sum(abs(temp).^2))
    
    %---------------------------------------------------------------------%
    % Plot the solutions                                                  %
    %---------------------------------------------------------------------%
   
     figure

     
    m_t = zeros(size(outTimes));

    errorTB = zeros(size(outTimes));
    
    for i = 1:length(outTimes)
       
        hold off
        
        subplot(1,2,1);
        
        rhot = Rho_t(i,:)';
        
        rhotTB = Rho_tTB(i,:)';
        
        rhotT = rhotTB(1:NT);
        rhotB = rhotTB(NT+1:NT+NB);
        
        rhoTI = InterpT.InterPol * rhotT;
        rhoBI = InterpB.InterPol * rhotB;
        
        rhoI = Interp.InterPol * rhot;
        
        tempBT = [reshape(rhoBI,NInterp2/2,NInterp1);reshape(rhoTI,NInterp2/2,NInterp1)];
        temp   = reshape(rhoI,NInterp2,NInterp1);
   
        errorTB(i) = sum(sum(abs(temp-tempBT).^2)) / sum(sum(abs(temp).^2));
        
        
        
        
        DCT.plot(rhotT);
       
        hold on
        
        DCB.plot(rhotB);
      
        hold on
        
        
        
        DC.plot(rhot);
             
        xlim([0;2]);
        ylim([0;1]);
        
        subplot(1,2,2);

        hold off
        
        plot(Pts.y1_kv(Ind.bottom),rhot(Ind.bottom),'b-o')
        
        hold on
        
        plot(PtsB.y1_kv(IndB.bottom),rhotB(IndB.bottom),'r--x')

        
        pause(0.005)
        
    end
    figure
    plot(outTimes, log10(errorTB))
    
%     for i=1:length(outTimes)
%         DC.plot(Rho_t(i,:)'); 
%         m_t(i) = Int*(Rho_t(i,:)');
%         title(['Interpolation of Solution at t = ', num2str(outTimes(i))]);
%         title(['mass = ', num2str(m_t(i))]);       
%         zlim([0 max(rho_ic)]);        
% 
%         pause(0.1);
% 
%     end   
%     
%     max(rho_ic - Rho_t(1,:)')
%     
%     figure
%     plot(outTimes,(m_t/m_ic -1 )*100,'o-')
    

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
        
        y = exp(-(y1-1).^2/0.1-(y2-0.25).^2)/0.1;
        %y = cos(y1*2*pi).*cos(y2*2*pi) + 1;
        
        %yb = exp(-(y1-1).^2-(y2-0.5).^2);
        %yt = exp(-(y1-1).^2-(y2-0.5).^2);
        %y = [yb; yt];
    end

end