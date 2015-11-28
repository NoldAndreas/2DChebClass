%Half Space


%Initialization
    vext  = @Vext15;        
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,...
                       'L1',2,'L2',2,'N',[20;20],...
                       'y2Min',0);

    PlotArea       = struct('y1Min',-5,'y1Max',5,'L1',3,'L2',2,...
                       'N2',100,'N1',100,'y2Min',0,'y2Max',4);    

    Phys_Area.Conv      = struct('L',3,'L2',1,'N',[30,30]);%'ep2Conv',0.1
    
    
    HS                             = HalfSpace(Phys_Area);%v2struct(L1,L2,N));
    HS.ComputeAll(PlotArea);     
    
    figure('color','white','Position',[0 0 250 250]);
    HS.PlotGridLines();    
    HS.PlotGrid();
    HS.PlotIsoline(0,'y2');
    HS.PlotIsoline(-1/sqrt(2),'y1');
    HS.PlotIsoline(1/sqrt(2),'y1');    
    ylim([-0.5 10]); xlim([-0.5 0.5]*10);
    	
    xlabel('$y_1$','Interpreter','Latex');
    ylabel('$y_2$','Interpreter','Latex');
    SaveFigure('HalfSpaceGrid');

%Composed HalfSpace

    N      = [20;20];
    N2bound = 10;
    L1 = 2;    L2 = 2;        
    y2Min = -0.5; h = 1;        
    vext  = @Vext15;       
    alpha = pi/2;%4;

    
    HSC         = HalfSpace_Composed(v2struct(L1,L2,N,N2bound,y2Min,h,alpha));            
    HSCCartPts  = HSC.GetCartPts();
    V           = vext(HSCCartPts.y1_kv,HSCCartPts.y2_kv);
    
    figure('color','white','Position',[0 0 250 250]);
    HSC.PlotGridLines(); 
    HSC.PlotGrid();
    
    HSC.Sub_HalfSpace.PlotIsoline(0,'y2');
    HSC.PlotIsoline(-1/sqrt(2),'y1');
    HSC.PlotIsoline(1/sqrt(2),'y1');	
    ylim([-0.5 10]); xlim([-0.5 0.5]*10);        	
    xlabel('$y_1$','Interpreter','Latex');
    ylabel('$y_2$','Interpreter','Latex');
    
    SaveFigure('HalfSpaceGrid_Composed');