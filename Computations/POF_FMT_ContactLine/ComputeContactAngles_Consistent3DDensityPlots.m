function ComputeContactAngles_Consistent3DDensityPlots()

   % Job_ComputeContactAngle(160,0.5,[-10 10],[0.5 20]);
   % Job_ComputeContactAngle(20,1.47,[0 30],[0.5 20]);
   bounds1    = [-10,10];
   bounds2    = [0.5,15];
   maxComp_y2 = 35;
   N          = [45,75];
   
   opts = v2struct(bounds1,bounds2,maxComp_y2,N);
   
   ComputAdditionalData([0.6,0.8,0.9,1.1,1.2,1.3,0.65,0.75,0.85,0.95,1.05,1.15,1.35]);
   
   
   opts.alpha_deg = 40;  opts.epw       = 1.375;   
   Job_ComputeContactAngle(opts);    
   
   opts.alpha_deg = 60;  opts.epw       = 1.25;   
   Job_ComputeContactAngle(opts);                
   
   opts.alpha_deg = 90;  opts.epw       = 1.0;      
   Job_ComputeContactAngle(opts);       
   
   opts.alpha_deg = 120;  opts.epw      = 0.7;   
   Job_ComputeContactAngle(opts); 
   
   opts.alpha_deg = 135;  opts.epw      = 0.55;   
   Job_ComputeContactAngle(opts); 
   
   %Job_ComputeContactAngle(20,1.47,bounds1,bounds2,25);    
   
    function ComputAdditionalData(epwList)
        for i = 1:length(epwList)
            opts.alpha_deg  = 90;  
            opts.epw        = epwList(i);
            opts.maxComp_y2 = 15;
            opts.N = [50,80];
            
            config = GetStandardConfig(opts);
            
            CLT = ContactLineHS(config);     
            CLT.Preprocess();
            CLT.ComputeEquilibrium();          	        
            %CLT.PostProcess();
            
            opts.alpha_deg = round(CLT.MeasureContactAngle(2,[10,14]));
            opts.N = [45,75];
            
            
            clear('CLT');
            close all;            
            
            
            config = GetStandardConfig(opts);
            
            CLT = ContactLineHS(config);     
            CLT.Preprocess();
            CLT.ComputeEquilibrium();          	        
            CLT.PlotDensityResult();
            
            clear('CLT');
            close all;            
            
        end
    end
end