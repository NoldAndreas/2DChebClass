function ThesisNanoscale_Fig5_6_7_8_9_ComputeContactAngles()

    
    AddPaths('ThesisNanoscale');   
    global dirData             
    
    %ComputeData(45,1.155,-2.5);
    %ComputeData(60,1.071,-5);
    %ComputeData(90,0.856,-7.5);
    %ComputeData(120,0.594,-7.5);
    %ComputeData(135,0.453,-10);
        
    PlotData(45,-2.5,'2015_9_19_19_55_33',[-0.1 0],[-0.15 0.05],0.6);
    PlotData(60,-5,'2015_9_20_0_2_34',[-0.2 -0.1 0],[-0.2 0.05],0.6);
    PlotData(90,-7.5,'2015_9_18_11_25_36',[-0.2 -0.1 0],[-0.2 0.05],0.55);
    PlotData(120,-7.5,'2015_9_20_4_1_22',[-0.1 0],[-0.1,0.05],0.48);
    PlotData(135,-10,'2015_9_20_7_53_56',[-0.04 -0.02 0],[-0.04 0],0.48);
    
    
    function ComputeData(alpha_deg,epw,bounds1)
        bounds1 = bounds1 + [0 15];
        AddPaths('ThesisNanoscale');   
        try
            opts = v2struct(alpha_deg,epw,bounds1);            
            opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);            
            Job_ComputeContactAngle(opts);
        catch err
            disp('ERROR')
            rethrow(err);        
        end

    end
    
    
    function PlotData(alpha_deg,bounds1,nameEq,yTicksDP,yLimsDP,subPlotPosX)        
        
        bounds1 = bounds1 + [0 15];
                
        AddPaths(['ThesisNanoscale' filesep 'deg' num2str(alpha_deg)]);            
        nameEq = [nameEq,'_'];

        
        try
                subplots = false;
                %**************************
                % ** Assemble plots **
                %**************************
                if(subplots)                
                    figMain = figure('Position',[0 0 500 600],'color','white');
                else
                    figMain = []; sP = [];
                    posSubPlotsGCA = [0.7 0.6 1.6 1.6];%250 200];
                    posSubPlots    = [0 0 2.4 2.5];%250 200];
                    units = 'inches';
                end

                %**************************
                % Upper left subplot: Contours
                if(subplots)                
                    sP = subplot(3,2,1);                                
                end                
                PutFigInSubplot([nameEq,'DensitySlices_contour'],figMain,sP);
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$y_2$','Interpreter','Latex');
                pbaspect([1 1 1]);
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(a)','LineWidth',0,'EdgeColor','none');
                
                if(~subplots)                                         
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);
                    SaveFigure(['SubFig1']);
                end

                % Upper right subplot: Density profiles
                if(subplots)                
                    sP = subplot(3,2,2);                                    
                end
                PutFigInSubplot([nameEq,'DensitySlices'],figMain,sP);
                xlim([0 15]);
                xlabel('$y_2$','Interpreter','Latex');
                ylabel('$\nDensity$','Interpreter','Latex');
                pbaspect([1 1 1]);
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(b)','LineWidth',0,'EdgeColor','none');
                if(~subplots)  
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);                    
                    SaveFigure(['SubFig2']);
                end

                % Middle left subplot: Contours
                if(subplots)                
                    sP = subplot(3,2,3);                    
                end
                PutFigInSubplot([nameEq,'DensitySlicesNormal_Contour'],figMain,sP);
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$y_2$','Interpreter','Latex');
                pbaspect([1 1 1]);
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(c)','LineWidth',0,'EdgeColor','none');
                if(~subplots)     
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);
                    SaveFigure(['SubFig3']);
                end

                % Middle right subplot: Density profiles
                if(subplots)                
                    sP = subplot(3,2,4);
                end
                PutFigInSubplot([nameEq,'DensityNormalInterface'],figMain,sP);
                xlim([-5 5]);
                xlabel('$y_2$','Interpreter','Latex');        
                ylabel('$\nDensity$','Interpreter','Latex');
                pbaspect([1 1 1]);
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(d)','LineWidth',0,'EdgeColor','none');
                if(~subplots)         
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);
                    SaveFigure(['SubFig4']);                    
                end

                % Bottom left subplot: Disjoining pressure
                if(subplots)                
                    sP = subplot(3,2,5);
                end
                PutFigInSubplot([nameEq,'DisjoiningPressures'],figMain,sP);                
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$\Pi$','Interpreter','Latex');
                pbaspect([1 1 1]);
                ylim(yLimsDP); set(gca,'YTick',yTicksDP);                
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(e)','LineWidth',0,'EdgeColor','none');
                if(~subplots)   
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);
                    SaveFigure(['SubFig5']);
                end

                % Bottom left subplot: Disjoining pressure
                close all;
                if(subplots)                
                    sP = subplot(3,2,6);                                    
                end
                figMain = PutFigInSubplot([nameEq,'AdsorptionIsotherm'],figMain,sP);
                %xlim(bounds1);
                xlabel('$\Delta \chemPot$','Interpreter','Latex');
                ylabel('$\ell$','Interpreter','Latex');
                ylim([0 15]);
                pbaspect([1 1 1]);                
                annotation('textbox',[0.02 0.9 0.1 0.05],'String','(f)','LineWidth',0,'EdgeColor','none');

                fsub  = openfig([dirData filesep nameEq,'FittingAdsorptionIsotherm']);
                xlabel('');  ylabel('');
                if(subplots)                         
                    inset2(figMain,fsub,0.1,[0.8,0.2]);               
                else
                    %inset2(figMain,fsub,0.5,[subPlotPosX,0.375]);               
                    inset2(figMain,fsub,0.4,[subPlotPosX,0.5]);               
                end
                close(fsub);                  
                
                if(subplots)
                    SaveFigure(['FullFigure_deg',num2str(alpha_deg)]);
                else
                    set(gcf,'Units',units,'Position',posSubPlots,'color','white');
                    set(gca,'Units',units,'Position',posSubPlotsGCA);
                    SaveFigure(['SubFig6']);    
                end
                close all;
          catch err
                disp('ERROR')
                rethrow(err);
                %msgString = getReport(exception,type,'hyperlinks',hlink);
                %disp(msgString);
            end
        
    end

    function f11 = PutFigInSubplot(filename,figMain,subP)
        f11  = openfig([dirData filesep filename]);
        if(~isempty(figMain))
            ax11 = gca;                                                

            f_cont = get(ax11,'children');

            set(0,'CurrentFigure',figMain)

            copyobj(f_cont,subP);
            close(f11);
        end
    end
    function config = GetStandardConfig(opts)                             
        
        config = ThesisNanoscale_GetStandardConfig(opts.alpha_deg,opts.epw);
                
        config.optsNum.PlotAreaCart       = struct('y1Min',opts.bounds1(1),'y1Max',opts.bounds1(2),...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'zMax',4,...
                                                    'N1',100,'N2',100);
    end
    function nameEq = Job_ComputeContactAngle(opts)

        config = GetStandardConfig(opts);
        close all;

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium(struct('solver','Picard'));      
        CLT.PostProcess(opts);
        CLT.PlotDensitySlices();
        CLT.PlotDensitySlicesNormalInterface();
        CLT.PlotDisjoiningPressures();                
        CLT.FittingAdsorptionIsotherm([10 14],1);
                
        [~,fn] = fileparts(CLT.FilenameEq);
        nameEq = [fn,'_'];            
    end
   
    function filename = ComputeExactAdsorptionIsotherm(opts)
        
        config = ThesisNanoscale_GetStandardConfig(90,opts.epw);
        
        config.optsNum.PhysArea.N         = [1,250];        
        config.optsNum.maxComp_y2         = -1;
        %config.optsNum.y1Shift            = 0;
        
        
        CL = ContactLineHS(config);
        CL.Preprocess();    close all;
        
        if(opts.alpha_deg > 90)
            optsDrying = 'drying';
        else
            optsDrying = 'wetting';
        end
        
        CL.ComputeAdsorptionIsotherm(750,optsDrying);
        
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(config.optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.3463);
        end
        
        filename = CL.AdsorptionIsotherm.Filename;

    end
end