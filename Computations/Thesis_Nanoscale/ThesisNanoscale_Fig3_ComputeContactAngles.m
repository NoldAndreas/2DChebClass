function ThesisNanoscale_Fig3_ComputeContactAngles()

    AddPaths('ThesisNanoscale');   
    global dirData 
        
    optsDefault = {'onlyPlot'};
    
    ComputeAndPlot(45,1.155,[0 20],optsDefault,'2015_9_19_19_55_33_');
    ComputeAndPlot(60,1.071,[-10 10],optsDefault);
    ComputeAndPlot(90,0.856,[-10 10],optsDefault);
    ComputeAndPlot(120,0.594,[-10 10],optsDefault);
    ComputeAndPlot(135,0.453,[-10 10],optsDefault);
    
    
    function ComputeAndPlot(alpha_deg,epw,bounds1,opts,nameEq)
        if(nargin < 4)
            opts = {};
        end
        
        if(IsOption(opts,'onlyPlot'))
            AddPaths(['ThesisNanoscale' filesep 'deg' num2str(alpha_deg)]);            
        else
            try
                opts = v2struct(alpha_deg,epw,bounds1);            
                opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
                nameEq = Job_ComputeContactAngle(opts);
            catch err
                disp('ERROR')
                rethrow(err);
                %msgString = getReport(exception,type,'hyperlinks',hlink);
                %disp(msgString);
            end
        end
        
        try
                %**************************
                % ** Assemble plots **
                %**************************
                figMain = figure('Position',[0 0 500 600],'color','white');

                %**************************
                % Upper left subplot: Contours
                s1 = subplot(3,2,1);
                PutFigInSubplot([nameEq,'DensitySlices_contour'],figMain,s1);
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$y_2$','Interpreter','Latex');
                pbaspect([1 1 1]);

                % Upper right subplot: Density profiles
                s2 = subplot(3,2,2);
                PutFigInSubplot([nameEq,'DensitySlices'],figMain,s2);
                xlim([0 15]);
                xlabel('$y_2$','Interpreter','Latex');
                ylabel('$\nDensity$','Interpreter','Latex');
                pbaspect([1 1 1]);

                % Middle left subplot: Contours
                s1 = subplot(3,2,3);
                PutFigInSubplot([nameEq,'DensitySlicesNormal_Contour'],figMain,s1);
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$y_2$','Interpreter','Latex');
                pbaspect([1 1 1]);

                % Middle right subplot: Density profiles
                s2 = subplot(3,2,4);
                PutFigInSubplot([nameEq,'DensityNormalInterface'],figMain,s2);
                xlim([-5 5]);
                xlabel('$y_2$','Interpreter','Latex');        
                ylabel('$\nDensity$','Interpreter','Latex');
                pbaspect([1 1 1]);

                % Bottom left subplot: Disjoining pressure
                s3 = subplot(3,2,5);
                PutFigInSubplot([nameEq,'DisjoiningPressures'],figMain,s3);
                xlim(bounds1);
                xlabel('$y_1$','Interpreter','Latex');
                ylabel('$\Pi$','Interpreter','Latex');
                pbaspect([1 1 1]);

                % Bottom left subplot: Disjoining pressure
                s4 = subplot(3,2,6);
                PutFigInSubplot([nameEq,'AdsorptionIsotherm'],figMain,s4);
                %xlim(bounds1);
                xlabel('$\Delta \chemPot$','Interpreter','Latex');
                ylabel('$\ell$','Interpreter','Latex');
                ylim([0 15]);
                pbaspect([1 1 1]);

                f11  = openfig([dirData filesep nameEq,'AdsorptionIsotherm']);
          catch err
                disp('ERROR')
                rethrow(err);
                %msgString = getReport(exception,type,'hyperlinks',hlink);
                %disp(msgString);
            end
        
    end

    function PutFigInSubplot(filename,figMain,subP)
        f11  = openfig([dirData filesep filename]);
        ax11 = gca;                                                
        
        f_cont = get(ax11,'children');
        
        set(0,'CurrentFigure',figMain)
                
        copyobj(f_cont,subP);
        close(f11);
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
        CL.ComputeAdsorptionIsotherm(750);
        
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(config.optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.3463);
        end
        
        filename = CL.AdsorptionIsotherm.Filename;

    end
end