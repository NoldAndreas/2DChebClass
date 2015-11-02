function ThesisNanoscale_Fig3_LiquidVapourDensityProfile

    AddPaths('ThesisNanoscale');       
    
    config = ThesisNanoscale_GetStandardConfig(90,1.0);    
    config.optsNum.PhysArea.N       = [60;3];
    config.optsNum.PhysArea.N2bound = 3;
	CL = ContactLineHS(config);
    CL.Preprocess();   
	
    %[~,pts1,rho_0_6] = GetLV(0.6);
    %[~,pts1,rho_0_7] = GetLV(0.7);
    %[~,pts1,rho_0_8] = GetLV(0.8);
    %[~,pts1,rho_0_9] = GetLV(0.9);    
    [~,pts1,rho_0_75] = GetLV(0.75);
    close all;
    figure('Position',[0 0 200 175],'color','white');
    %plot(pts1,rho_0_6); hold on;
    %plot(pts1,rho_0_7);
    %plot(pts1,rho_0_8);
    %plot(pts1,rho_0_9);
    plot(pts1,rho_0_75,'k'); hold on;
    plot([-5 5],CL.optsPhys.rhoLiq_sat*[1 1],'k--');
    plot([-5 5],CL.optsPhys.rhoGas_sat*[1 1],'k--');
    
    xlabel('$y_1$','Interpreter','Latex');
    ylabel('$\nDensity$','Interpreter','Latex');
    xlim([-5 5]);
    ylim([0 0.8]);
	pbaspect([1 1 1]);
    SaveFigure(['LiquidVapourInterface']);
        

    function [rho1D,pts_IP,rho_IP] = GetLV(T)
        CL.SetkBT(T);        
        
        close all;
        [om,rho1D,params] = CL.Compute1D('LV');            
        
        opts = struct('dist0',true,'plain',true);
                
        [pts_IP,rho_IP] = CL.IDC.do1DPlotParallel(rho1D); close all;        
        %title(['T = ',num2str(T)],'Interpreter','Latex');        
        %ylabel('$\nDensity$','Interpreter','Latex');        
        %plot(CLN.IDC.Pts.y2,rho1D); 
        %xlim([0 10]);
    end
end