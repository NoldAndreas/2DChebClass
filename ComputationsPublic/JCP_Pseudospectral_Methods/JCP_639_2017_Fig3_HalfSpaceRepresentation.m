function JCP_639_2017_Fig3_HalfSpaceRepresentation()

    AddPaths('JCP_639_2017');
    close all;

    %*************************************
    % (1) Representation of half space
    %*************************************

    % **** Initialization *****
    Phys_Area       = struct('N',[20;20],...
                             'L1',2,'L2',2,...
                             'y2Min',0);
    HS              = HalfSpace(Phys_Area);     
    
    % **** Plotting ****
    figure('color','white','Position',[0 0 600 600]);
    HS.PlotGridLines();    
    HS.PlotGrid();
    HS.PlotIsoline(0,'y2');
    HS.PlotIsoline(-1/sqrt(2),'y1');
    HS.PlotIsoline(1/sqrt(2),'y1');    
    ylim([-0.5 10]); xlim([-5 5]);    	
    xlabel('$y_1$','Interpreter','Latex','fontsize',25);
    ylabel('$y_2$','Interpreter','Latex','fontsize',25);
    
    SaveFigure('JCP_639_2017_Fig3_a');


    %********************************************
    % (2) Representation of composed half space
    %********************************************

    Phys_Area   = struct('N',[20;20],'N2bound',10,...        
                         'L1',2,'L2',2,...
                         'y2Min',-0.5,'h',1);    
    HSC         = HalfSpace_Composed(Phys_Area);            
    
    % **** Plotting ****
    figure('color','white','Position',[0 0 600 600]);
    HSC.PlotGridLines(); 
    HSC.PlotGrid();    
    HSC.Sub_HalfSpace.PlotIsoline(0,'y2');
    HSC.PlotIsoline(-1/sqrt(2),'y1');
    HSC.PlotIsoline(1/sqrt(2),'y1');
	ylim([-0.5 10]); xlim([-5 5]);
    xlabel('$y_1$','Interpreter','Latex','fontsize',25);
    ylabel('$y_2$','Interpreter','Latex','fontsize',25);
    
    SaveFigure('JCP_639_2017_Fig3_HalfSpaceRepresentation_b');
    
end