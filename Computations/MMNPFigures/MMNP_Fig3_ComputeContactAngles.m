function MMNP_Fig3_ComputeContactAngles()

   bounds1    = [-10,10];
   bounds2    = [0.5,15];
   maxComp_y2 = 35;
   N          = [45,75];
   
   opts = v2struct(bounds1,bounds2,maxComp_y2,N);
      
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
   


   % Job_ComputeContactAngle(160,0.5,[-10 10],[0.5 20]);
   % Job_ComputeContactAngle(20,1.47,[0 30],[0.5 20]);
   %Job_ComputeContactAngle(21,1.47,[-5 30],[0.5 15],25);    
   
%    Job_ComputeContactAngle(90,1.0,[-10 10],[0.5 15]); 
%    Job_ComputeContactAngle(135,0.55,[-10 10],[0.5 15]); 
%    Job_ComputeContactAngle(40,1.375,[0 20],[0.5 18]);
%    Job_ComputeContactAngle(120,0.7,[-10 10],[0.5 18]);      
%    Job_ComputeContactAngle(60,1.25,[-10 10],[0.5 18]);   
    
end