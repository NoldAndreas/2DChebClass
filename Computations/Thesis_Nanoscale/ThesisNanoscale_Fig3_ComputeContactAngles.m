function ThesisNanoscale_Fig3_ComputeContactAngles()

   AddPaths('ThesisNanoscale');   
   

    try 
        opts.bounds1   = [-10,10];
       opts.alpha_deg = 90;  
       opts.epw       = 0.856;     
       opts.dryingWetting = 'wetting';
       opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
       Job_ComputeContactAngle(opts); 
    catch err
        disp('ERROR')
    end

    try            
       opts.bounds1   = [0,20];
       opts.alpha_deg = 45;  
       opts.epw       = 1.155;   
       opts.dryingWetting = 'wetting';   
       opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
       Job_ComputeContactAngle(opts);    
    catch err
        disp('ERROR')
    end       
       
    try
       opts.bounds1   = [-10,10];
       opts.alpha_deg = 60;  
       opts.epw       = 1.071;
       opts.dryingWetting = 'wetting';
       opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);   
       Job_ComputeContactAngle(opts);                      
    catch err
        disp('ERROR')
    end       
   
    try
       opts.alpha_deg = 120;  
       opts.epw       = 0.594;
       opts.dryingWetting = 'drying';
       opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
       Job_ComputeContactAngle(opts); 
	catch err
        disp('ERROR')
    end
   
    try
       opts.alpha_deg = 135;  
       opts.epw       = 0.453;
       opts.dryingWetting = 'drying';
       opts.AdsorptionIsotherm_file = ComputeExactAdsorptionIsotherm(opts);
       Job_ComputeContactAngle(opts);
	catch err
        disp('ERROR')
    end

    function Job_ComputeContactAngle(opts)

        config = GetStandardConfig(opts);
        close all;

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium(struct('solver','Picard'));      
        CLT.PostProcess(opts);
        CLT.PlotDensitySlices();
        CLT.PlotDisjoiningPressures();                
    end
    function config = GetStandardConfig(opts)                             
        
        config = ThesisNanoscale_GetStandardConfig(opts.alpha_deg,opts.epw);
                
        config.optsNum.PlotAreaCart       = struct('y1Min',opts.bounds1(1),'y1Max',opts.bounds1(2),...
                                                   'y2Min',0.5,'y2Max',15.5,...
                                                   'zMax',4,...
                                                    'N1',100,'N2',100);
    end
    function filename = ComputeExactAdsorptionIsotherm(opts)
        
        config = ThesisNanoscale_GetStandardConfig(90,opts.epw);
        
        config.optsNum.PhysArea.N         = [1,250];        
        config.optsNum.maxComp_y2         = -1;
        %config.optsNum.y1Shift            = 0;
        
        
        CL = ContactLineHS(config);
        CL.Preprocess();    close all;
        CL.ComputeAdsorptionIsotherm(500,opts.dryingWetting);    
        
        CL.FittingAdsorptionIsotherm([10 14],1)
        if(config.optsPhys.kBT == 0.75)
            CL.SumRule_AdsorptionIsotherm(0.3463);
        end
        
        filename = CL.AdsorptionIsotherm.Filename;

    end
end