function ThesisNanoscale_LiquidDensityProfile

    AddPaths('ThesisNanoscale');       
    
    config = ThesisNanoscale_GetStandardConfig(90,1.0);    
    config.optsNum.PhysArea.N = [1;100];
	CLN = ContactLineHS(config);
    CLN.Preprocess();   
	
    GetWL(0.6);
    GetWL(0.7);
    GetWL(0.8);
    GetWL(0.9);
    close all;

    function [rho1D] = GetWL(T)
        CLN.SetkBT(T);        
        
        close all;
        [om,rho1D,params] = CLN.Compute1D('WL');            
        
        opts = struct('dist0',true,'plain',true);
        
        figure('Position',[0 0 175 135],'color','white');
        CLN.IDC.plotLine([0 0],[0.5 5.5],rho1D,opts);   
        title(['T = ',num2str(T)],'Interpreter','Latex');
        xlabel('$y_2$','Interpreter','Latex');
        %ylabel('$\nDensity$','Interpreter','Latex');        
        xlim([0 5]);
        ylim([0 4]);
        pbaspect([1 1 1]);
        %annotation('textbox',[0.25 0.85 0.4 0.05],'String',['T = ',num2str(T)],'LineWidth',0,'EdgeColor','none');
        SaveFigure(['WallLiquidInterface_T_',num2str(T)]);
        %plot(CLN.IDC.Pts.y2,rho1D); 
        %xlim([0 10]);
    end
end