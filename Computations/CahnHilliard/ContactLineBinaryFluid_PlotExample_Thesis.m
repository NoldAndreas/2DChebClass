function ContactLineBinaryFluid_PlotExample_Thesis
close all;
    AddPaths('DIBinaryPaper');            

    PhysArea = struct('N',[50,40],'y2Min',0,'y2Max',66,...
                      'L1',7,'IntInterval',[-10,10]);%,'NBorder',[30,200,30,200]);

    PlotArea = struct('y1Min',-15,'y1Max',15,'N1',80,'N2',80);   
    SubArea  = struct('shape','Box','N',[60,60],...
                      'y1Min',-5,'y1Max',5,'y2Min',0,'y2Max',10);   	
                  
    %optsSolv = struct('noIterations',40,'lambda',0.8,'Seppecher_red',1);

    optsNum  = v2struct(PhysArea,PlotArea);                   	

    optsPhys = struct('thetaEq',pi/2,...       
                       'theta',90*pi/180,...
                       'Cak',0.005,'Cn',1,...
                       'UWall',1,...
                       'l_diff',2.75,...%'mobility',10,...
                       'nParticles',0);

    config = v2struct(optsPhys,optsNum);    
        
    %******************************************
    %******************************************
    %******************************************
        
    Cn                            = 1/config.optsPhys.l_diff;  
    Cak                           = config.optsPhys.Cak;
    

    xL = [-7 7]; yL = [0 14];
    config.optsNum.PlotArea = struct('y1Min',xL(1)/Cn,'y1Max',xL(2)/Cn,...
                                     'y2Min',yL(1)/Cn,'y2Max',yL(2)/Cn,'N1',80,'N2',80);   

    DI = DiffuseInterfaceBinaryFluid(config);
    DI.Preprocess();
    DI.IterationStepFullProblem();                    
    DI.PostProcess();        

    PtsCart = DI.IDC.GetCartPts;

    %Plot Values at stagnation point    
    L    = 15;
    st0  = GetValuesParallelToWall(0,L); hold on;
    stSP = GetValuesParallelToWall(DI.StagnationPoint.y2_kv(1),L);
        
    figure('Position',[0 0 400 400],'color','white');
    PlotValuesParallelToWall(st0,'-',DI.IsoInterface.h(1));
    PlotValuesParallelToWall(stSP,'--',interp1(DI.IsoInterface.y2,DI.IsoInterface.h,DI.StagnationPoint.y2_kv(1)));%DI.StagnationPoint.y1_kv(1));

    xlim([-L,L]);
    set(gca,'linewidth',1.5); set(gca,'fontsize',20);
    xlabel('$y_1$','fontsize',20,'Interpreter','Latex');       
    pbaspect([1 1 1]);        
    config.stagnationPoint = DI.StagnationPoint;
	SaveFigure('ValuesThroughStagnationPoint',config);               

%     y1    = [-10,10]/Cn+DI.StagnationPoint.y1_kv(1);
%     [mu,pts] = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.mu); close all; mu = mu/(Cn*Cak);
%     [u,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1:end/2)); close all;
%     [v,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1+end/2:end)); close all;
%     [p,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.p); close all; p = p/Cn;
%     [j1,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1:end/2)); close all;
%     [j2,pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1+end/2:end)); close all;
% 
%     lw  = 1.5;
%     y   = (pts.y1_kv - DI.StagnationPoint.y1_kv(1))*Cn;
%     figure('Position',[0 0 800 800],'color','white');
%     plot(y,mu,'-k','linewidth',lw); hold on;
%     plot(y,p,'--k','linewidth',lw);
%     plot(y,u,'-b','linewidth',lw);
%     plot(y,v,':b','linewidth',lw);
%     plot(y,j1,'-m','linewidth',lw);
%     plot(y,j2,':m','linewidth',lw);
% 
%     set(gca,'linewidth',1.5); set(gca,'fontsize',20);
%     xlabel('$y_1$','fontsize',20,'Interpreter','Latex');       
%     pbaspect([1 1 1]);        
%    config.stagnationPoint = DI.StagnationPoint;

    %Plot velocities
    figure('Position',[0 0 800 800],'color','white');       
    DI.PlotResultsPhi();
    DI.PlotU();
    DI.PlotStagnationPoint();
   % DI.IDC.plotFlux(DI.uv,[],[],1,'k');%[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = 
   % DI.IDC.plotFlux(DI.flux); %[y1_s,y2_s,fl_y1_q,fl_y2_q] =         
    SetTicksLabels(Cn,[-5,0,5],[0,5,10]);                 
    SaveFigure('StagnationPoint_Velocity',config);

    %Plot flux j_1-j_2 next to stagnation point
    figure('Position',[0 0 800 800],'color','white');       
    DI.PlotResultsPhi();
    DI.PlotU(DI.flux); %[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = 
    DI.PlotStagnationPoint();
    %DI.IDC.plotFlux(DI.flux); %[y1_s,y2_s,fl_y1_q,fl_y2_q] =         
    SetTicksLabels(Cn,[-5,0,5],[0,5,10]);                                 
    SaveFigure('StagnationPoint_Flux',config);

    %******* Plot 3D chemical potential, pressure and phi
    vals = {'mu','p','phi'};
    labs = {'$\mu$','p','$\phi$'};
    for k = 1:3            
        figure('Position',[0 0 800 800],'color','white');        
        if(strcmp(vals{k},'mu'))
             DI.IDC.plot(DI.(vals{k})/(Cn*Cak)  );
        elseif(strcmp(vals{k},'p'))
             DI.IDC.plot(DI.(vals{k})/Cn  );
        else
             DI.IDC.plot(DI.(vals{k}));                 
        end                        
        DI.PlotU([],struct('color','b','linewidth',1.4));%close all;[y1M,y2M,VIM] =   %close all;[y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] =            
        SetTicksLabels(Cn,[-5,0,5],[0,5,10]);
       % surf(y1M*Cn,y2M*Cn,VIM);  hold on;
       % h = streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   
       % xlim(xL); ylim(yL);
       % pbaspect([(xL(2)-xL(1))/(yL(2)-yL(1)) 1 1]);
       % set(h,'linewidth',1.5);  %set(h,'color','k');
        %set(gca,'linewidth',1.5); set(gca,'fontsize',20);
        %xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
        %ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
        zlabel(labs{k},'fontsize',25,'Interpreter','Latex');  
        SaveFigure([vals{k},'_3DPlot'],config);
    end  


     %         %Plot velocities next to stagnation point
%         [y1M,y2M,VIM] = DI.IDC.plot(DI.phi,'contour'); close all;
%         [y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = DI.PlotU(); close all;                
%         [y1_s,y2_s,fl_y1_q,fl_y2_q] = DI.IDC.plotFlux(DI.uv); close all;
%         
% 
%         figure('Position',[0 0 800 800],'color','white');
%         [C,h] = contour(y1M*Cn,y2M*Cn,VIM,'linewidth',2.5);  hold on;        
%         h = streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   hold on;
%         quiver(y1_s*Cn,y2_s*Cn,fl_y1_q,fl_y2_q);        
%         %quiver(y1_s_Sepp,y2_s_Sepp,fl_y1_q__Sepp,fl_y2_q__Sepp,'color','m');        
%         plot(DI.IsoInterface.h*Cn,DI.IDC.Pts.y2*Cn,'k','linewidth',3);
%         plot(DI.StagnationPoint.y1_kv(1)*Cn,DI.StagnationPoint.y2_kv(1)*Cn,'o','MarkerFaceColor','m','MarkerSize',12)
%         
%         set(gca,'linewidth',1.5); set(gca,'fontsize',20);
%         xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
%         ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
%         pbaspect([1 1 1]);
%         SaveFigure('StagnationPoint',config);

    %************************************
    %************************************
    %************************************
    %Plot wider area
    DI.optsNum.PlotArea = struct('y1Min',-10/Cn,'y1Max',10/Cn,...
                                 'y2Min',0/Cn,'y2Max',20/Cn,'N1',80,'N2',80);   
    DI.IDC.InterpolationPlotCart(DI.optsNum.PlotArea,true);

     %Plot flux j_1-j_2 next to stagnation point
     %Plot velocities next to stagnation point         
     uSepp = GetSeppecherSolutionCart_Blurred([PtsCart.y1_kv - DI.deltaX,...
                                            PtsCart.y2_kv],1,0,0,DI.theta);
    [y1_s_Sepp,y2_s_Sepp,fl_y1_q__Sepp,fl_y2_q__Sepp] = DI.IDC.plotFlux(uSepp); close all;                                                             
    [y1MU_Sepp,y2MU_Sepp,fl_y1_Sepp,fl_y2_Sepp,startMask1_Sepp,startMask2_Sepp] = DI.PlotU(uSepp); close all;                

    [y1M,y2M,VIM] = DI.IDC.plot(DI.phi,'contour'); close all;
    [y1MU,y2MU,fl_y1,fl_y2,startMask1,startMask2] = DI.PlotU(); close all;                
    [y1_s,y2_s,fl_y1_q,fl_y2_q] = DI.IDC.plotFlux(DI.uv); close all;


    figure('Position',[0 0 800 800],'color','white');
    [C,h] = contour(y1M*Cn,y2M*Cn,VIM,'linewidth',2.5);  hold on;        
    streamline(y1MU*Cn,y2MU*Cn,fl_y1,fl_y2,startMask1*Cn,startMask2*Cn);   hold on;
    h = streamline(y1MU_Sepp*Cn,y2MU_Sepp*Cn,fl_y1_Sepp,fl_y2_Sepp,startMask1_Sepp*Cn,startMask2_Sepp*Cn);   hold on;
    set(h,'color','m');
    quiver(y1_s*Cn,y2_s*Cn,fl_y1_q,fl_y2_q);        
    quiver(y1_s_Sepp*Cn,y2_s_Sepp*Cn,fl_y1_q__Sepp,fl_y2_q__Sepp,'color','m');        
    plot(DI.IsoInterface.h*Cn,DI.IDC.Pts.y2*Cn,'k','linewidth',3);
    plot(DI.StagnationPoint.y1_kv(1)*Cn,DI.StagnationPoint.y2_kv(1)*Cn,'o','MarkerFaceColor','m','MarkerSize',12)

    set(gca,'linewidth',1.5); set(gca,'fontsize',20);
    xlabel('$y_1$','fontsize',20,'Interpreter','Latex');
    ylabel('$y_2$','fontsize',20,'Interpreter','Latex');
    pbaspect([1 1 1]);
    SaveFigure('CoxComparison',config);                  
    
    function st = GetValuesParallelToWall(y2_SP,L)

        %Plot Values at stagnation point
        y1       = [-L,L]/Cn+DI.StagnationPoint.y1_kv(1);
        [st.phi] = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.phi); close all; 
        [mus]    = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.mu); close all; 
        st.mu    = mus/(Cn*Cak);
        [st.u]   = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1:end/2)); close all;
        [st.v]   = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.uv(1+end/2:end)); close all;
        [p]      = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.p); close all; 
        st.p     = p/Cn;
        [st.j1]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1:end/2)); close all;
        [st.j2,st.pts]  = DI.IDC.plotLine(y1,y2_SP*[1 1],DI.flux(1+end/2:end)); close all;

    end
    function PlotValuesParallelToWall(st,style,y1Ref)
        lw  = 1.5;
        y   = (st.pts.y1_kv - y1Ref)*Cn;        
        plot(y,st.mu,style,'linewidth',lw,'color','m'); hold on;
        plot(y,st.p,style,'linewidth',lw,'color','g');
        plot(y,st.u,style,'linewidth',lw,'color','b');
        plot(y,st.v,style,'linewidth',lw,'color','c');
        plot(y,st.phi,style,'linewidth',lw,'color','k');
        %plot(y,st.j1,'-m','linewidth',lw);
        %plot(y,st.j2,':m','linewidth',lw);
    end

	function SetTicksLabels(Cn,xlRed,ylRed)
        
        xl    = xlim;
        if(nargin < 2)
            xlRed = xl(1)  + [0.1,0.5,0.9]*(xl(2)-xl(1));
            xlRed = xlRed*Cn;
        end
        for i = 1:length(xlRed)
            xlRed(i) = round(xlRed(i))/Cn;
        end
        set(gca,'XTick',xlRed);                               
        
        T = get(gca,'Xtick');                
        c = {};
        R = 100;
        
        for i = 1:length(T)
            c{end+1} = num2str(round(R*T(i)*Cn)/R);
        end
        set(gca,'XTickLabel',c);
        
        
        yl    = ylim;
        if(nargin < 3)
            ylRed = yl(1)  + [0,0.5,0.9]*(yl(2)-yl(1));
            ylRed    = ylRed*Cn;
        end        
        for i = 1:length(ylRed)
            ylRed(i) = round(ylRed(i))/Cn;
        end
        set(gca,'YTick',ylRed);        
        
        T = get(gca,'Ytick');                
        c = {};
        for i = 1:length(T)
            c{end+1} = num2str(round(R*T(i)*Cn)/R);
        end
        set(gca,'YTickLabel',c);                          
        
        set(gca,'linewidth',1.5); set(gca,'fontsize',25);
        xlabel('$y_1$','fontsize',25,'Interpreter','Latex');
        ylabel('$y_2$','fontsize',25,'Interpreter','Latex');        
        pbaspect([1 1 1]);
     
    end

end