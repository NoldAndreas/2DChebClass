function ComputeContactAngles_Consistent3DDensityPlots()

   % Job_ComputeContactAngle(160,0.5,[-10 10],[0.5 20]);
   % Job_ComputeContactAngle(20,1.47,[0 30],[0.5 20]);
   bounds1    = [-10,10];
   bounds2    = [0.5,15];
   maxComp_y2 = 35;
   N          = [45,75];
   
   opts    = v2struct(bounds1,bounds2,maxComp_y2,N);
   
    
   
   epwList   = [0.8,0.6,0.9,1.1,1.2,1.3,0.65,0.75,0.85,0.95,1.05,1.15,1.35];   
   
   epwList1  = [1.375,1.35,1.32,1.285,1.25,1.211,1.1715,1.13,1.1,0.7,0.65,0.55];%,1.407,1.431,];   
   alpha1    = [40,45,50,55,60,65,71,76,79,120,125,135];
   
   epwList3 = [1.1715,1.407,1.431,];
   
    for i = 9%:length(epwList1)
        ComputeEpw(alpha1(i),epwList1(i));   
      %  alpha_deg = GetAlphaDeg(epwList1(i));
      %  ComputeEpw(alpha_deg,epwList1(i));      
    end
   ComputeEpw(alpha1(10),epwList1(10));   
  
   
   %epwList2  = [0.75,0.8,0.85,0.9];
%    alphaDeg2 = [115,110,105,100];
%   
%    for i = 2:length(epwList2)
%         ComputeEpw(alphaDeg2(i),epwList2(i));      
    %end

   for i = 1:length(epwList3)
        alpha_deg = GetAlphaDeg(epwList3(i));
        ComputeEpw(alpha_deg,epwList3(i));      
   end
   for i = 1:length(epwList1)
        ComputeEpw(alpha1(i),epwList1(i));   
      %  alpha_deg = GetAlphaDeg(epwList1(i));
      %  ComputeEpw(alpha_deg,epwList1(i));      
   end
   for i = 1:length(epwList)
        alpha_deg = GetAlphaDeg(epwList(i));
        ComputeEpw(alpha_deg,epwList(i));      
   end
   
   
      
 %  ComputeEpw(40,1.375);       
 %  ComputeEpw(60,1.25);       
 %  ComputeEpw(90,1.0);
 %  ComputeEpw(120,0.7);
 %  ComputeEpw(135,0.55);
   
%    opts.alpha_deg = 60;  opts.epw       = 1.25;   
%    Job_ComputeContactAngle(opts);                
%    
%    opts.alpha_deg = 90;  opts.epw       = 1.0;      
%    Job_ComputeContactAngle(opts);       
%    
%    opts.alpha_deg = 120;  opts.epw      = 0.7;   
%    Job_ComputeContactAngle(opts); 
%    
%    opts.alpha_deg = 135;  opts.epw      = 0.55;   
%    Job_ComputeContactAngle(opts); 
   
   %Job_ComputeContactAngle(20,1.47,bounds1,bounds2,25);    
   
    function ComputeEpw(alpha_deg,epw)
        opts.alpha_deg  = alpha_deg;  
        opts.epw        = epw;
        opts.N          = [45,75];
        opts.maxComp_y2 = 15;        

        config = GetStandardConfig(opts);

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium();          	        
        CLT.PostProcess();    
                    
        CLT.PlotDensityResult('DP'); 
        
        clear('CLT');
        close all;    
    end
    function alpha_deg = GetAlphaDeg(epw)
        opts.alpha_deg  = 90;  
        opts.epw        = epw;
        opts.maxComp_y2 = 15;
        opts.N = [50,80];

        config = GetStandardConfig(opts);

        CLT = ContactLineHS(config);     
        CLT.Preprocess();
        CLT.ComputeEquilibrium();          	        
        %CLT.PostProcess();        

        alpha_deg = round(CLT.MeasureContactAngle(2,[10,14]));
        close all;
        clear('CLT');
    end
end