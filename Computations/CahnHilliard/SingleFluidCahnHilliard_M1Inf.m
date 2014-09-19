 function SingleFluidCahnHilliard_M1Inf()

    close all;
    
    global dirData
    AddPaths();        
    ChangeDirData([dirData filesep 'CahnHilliard_InnerRegion'],'ORG');    
    %% Parameters    
    %PhysArea = struct('N',[70,40],'y2Min',0,'y2Max',20,'L1',7,... %12,80,50
    PhysArea = struct('N',[50,30],'y2Min',0,'y2Max',10,'L1',7,... %12,80,50    
                      'NBorder',200,...
                      'y1Max',inf);

    PlotArea = struct('y1Min',-20,'y1Max',20,'N1',80,...
                      'y2Min',0,'y2Max',PhysArea.y2Max,'N2',80);   

	optsNum  = v2struct(PhysArea,PlotArea);                   	
    
    optsPhys = struct('thetaEq',pi/2,...   
                      'zeta',10,'eta',1,...
                      'Cak',0.1,'Cn',1,...
                      'UWall',1,...
                      'phi_m',4,...
                      'nParticles',0,...
                      'fluidInterface','top');
                    
    config = v2struct(optsPhys,optsNum);   
    
    y2M = [(10:2.5:15)';(17.5:2.5:27.5)'];
    stagnationPoint = zeros(length(y2M),2);
                      
    for j = 1:length(y2M)
        config.optsNum.PhysArea.y2Max = y2M(j);
        config.optsNum.PlotArea.y2Max = y2M(j);
        DI = DiffuseInterfaceSingleFluid(config);
        DI.Preprocess();        
        if(y2M(j)>15)
            DI.SolveMovingContactLine(30);    
        else
            DI.SolveMovingContactLine(25);    
        end
        DI.SavePlotResults();
        DI.PlotErrorIterations();
        
        stagnationPoint(j,:) = DI.StagnationPoint;
        clear 'DI'
        close all;
    end
    
    figure('color','white','Position',[0 0 800 800]);
    plot(y2M,stagnationPoint(:,2),'o','MarkerFaceColor','k','MarkerSize',10); hold on;
    plot(y2M,stagnationPoint(:,1),'d','MarkerFaceColor','k','MarkerSize',10);

    plot(y2M,stagnationPoint(1,2) + (y2M-y2M(1))*(stagnationPoint(end,2)-stagnationPoint(1,2))/(y2M(end)-y2M(1)),'r','linewidth',2);
    plot(y2M,stagnationPoint(1,1) + (y2M-y2M(1))*(stagnationPoint(end,1)-stagnationPoint(1,1))/(y2M(end)-y2M(1)),'r','linewidth',2);

    xlabel('$y_{2,max}$','Interpreter','Latex','fontsize',20);
    ylabel('$y_{1/2,SP}$','Interpreter','Latex','fontsize',20);
    set(gca,'fontsize',17.5);
    set(gca,'linewidth',1.5);
    
	print2eps([dirData filesep 'StagnationPoints' ],gcf);
	saveas(gcf,[dirData filesep 'StagnationPoints.fig']);
    


end
