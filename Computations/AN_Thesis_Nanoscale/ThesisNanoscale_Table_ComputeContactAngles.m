function ThesisNanoscale_Table_ComputeContactAngles()

    
    AddPaths('ThesisNanoscale');               
    
    ComputeData(45,1.155,-2.5);
    ComputeData(60,1.071,-5);
    ComputeData(90,0.856,-7.5);
    ComputeData(120,0.594,-7.5);
    ComputeData(135,0.453,-10);
        
    
    function ComputeData(alpha_deg,epw,bounds1)
        %bounds1 = bounds1 + [0 15];
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
        
    function config = GetStandardConfig(opts)                             
        
        config = ThesisNanoscale_GetStandardConfig(opts.alpha_deg,opts.epw);
                
        config.optsNum.PlotAreaCart       = struct('y1Min',opts.bounds1(1),'y1Max',opts.bounds1(2),...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'zMax',4,...
                                                    'N1',100,'N2',100);
    end
    function st = Job_ComputeContactAngle(opts)

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
        
        st.epw     = 
        st.thYoung =         
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