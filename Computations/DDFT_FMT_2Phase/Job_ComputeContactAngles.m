function Job_ComputeContactAngles()

    global dirData
    AddPaths();    
    
    ChangeDirData([dirData filesep 'FMT_CLEq_BH_40X40_epw'],'ORG');    

    PhysArea = struct('N',[40,65],...
                      'L1_Skewed',5,...
                      'L2_Skewed',2,...
                      'L2_AD_Skewed',2.,...
                      'y2wall',0.,...
                      'N2bound',24,'h',1,...
                      'alpha_deg',14);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[34,34]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',50,'N2disc',50);                   
%   Fex_Num   = struct('Fex','CarnahanStarling');
 
    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',20,...
                     'y1Shift',0);

    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.49);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                        
    
    %***********************************************************
    %Setup result file for this Job
    filename   = ([dirData filesep]); 
    filename   = [filename,'Job__12_11_13_ComputeContactAngles_epw_',...
                                getTimeStr(),'.txt'];
    
    %filename    = [dirData filesep subDir filesep 'Job__12_11_13_ComputeContactAngles_epw.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
    
    %epw    = [(0.74:0.05:1.09),1.1,1.11,1.12];
    %degAngle    = 90;
    
 %   epw       = (1.46:0.02:1.6);
 %   degAngle  = 30.;
    
 %   epw       = (1.48:0.01:1.55);
 %   degAngle  = 19.;
    
    epw = 1.49;  degAngle = 14; %FMT
    epw = 2.58;  degAngle = 14; %LDA
    epw = 2.57;  degAngle = 17; %LDA
    
    epw = 1.49;
    %LDA: 
    %epsilonw - angle
    %1.5      - 88
    %2        - 60
    %2.4      - 33
    %2.57     - 15.5

    config.optsNum.PhysArea.alpha_deg = degAngle;      
    
    
    for i = 1:length(epw)
        config.optsPhys.V1.epsilon_w  = epw(i);
        %**************************************************************************
        %Step 1: Test accuracy of computing surface tensions and deviation 
        %         in contact angle computation
        %**************************************************************************            
        ConvergenceSurfaceTensions(config);
        
        close all;
        ChangeDirData();                
                                
         CL = ContactLine(config); %CL with normal high resolution
         %CL.Preprocess();        
         %CL.ComputeST();
         CL.Compute();
%         
%         %Print expected error in contact angle
%         alpha_HIRES = ComputeContactAngle(CL_N.ST_1D.om_wallGas,...
%                                           CL_N.ST_1D.om_wallLiq,...
%                                           CL_P.ST_1D.om_LiqGas);
%         fprintf(['High resolution CA: ',num2str(alpha_HIRES*180/pi,3),'[deg] ,',...
%                 'currect resolution CA: ',num2str(CL.alpha_YCA*180/pi,3),'[deg] ,',...
%                 'Error: ',num2str((alpha_HIRES-CL.alpha_YCA)*180/pi,3),' [deg].\n']);
%         
%         
%         sol = ThreePhaseContactLine_FMT_BH(config);
%         
%         config.optsNum.PhysArea.N = n;
       %**************************************************************************
       %**************************************************************************        
%         fileID = fopen(filename,'a');
%         fprintf(fileID,'# ******************* \n');
%         fprintf(fileID,'# %s \n',['epw =  ',num2str(epw(i),3)]);
%         fclose(fileID);
%         
%          %(0) Compute Young CA for a 90 [deg] grid
%          config.optsNum.maxComp_y2         = -1; %No 2D Computation
%          config.optsNum.PhysArea.alpha_deg = 90; 
%          sol = ThreePhaseContactLine_FMT_BH(config);
%          fileID      = fopen(filename,'a');
%          fprintf(fileID,'# %s \n',['Young CA for 90 deg grid: ',num2str(sol.Young_CA*180/pi,3)]);
%          fclose(fileID);
%         
%         
%         %(1) 1st Iteration
%        % config.optsNum.PhysArea.alpha_deg = round(degAngle);
%         config.optsNum.maxComp_y2         = 12;
%         
%         sol                               = ThreePhaseContactLine_FMT_BH(config,subDir);                        
%         
%         fileID      = fopen(filename,'a');
%         fprintf(fileID,'# %s \n','** 1st Iteration **');
%         fprintf(fileID,'# %s \n',['maxComp_y2: ',num2str(config.optsNum.maxComp_y2,3)]);
%         fprintf(fileID,'# %s \n',['CA Grid: ',num2str(config.optsNum.PhysArea.alpha_deg,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
%         fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
%         fclose(fileID);
% %         
%         %(2) 2nd Iteration
%         config.optsNum.PhysArea.alpha_deg = round(sol.Measured_StaticCA*180/pi); %update grid CA
%         config.optsNum.maxComp_y2         = 16;
%         
%         sol = ThreePhaseContactLine_FMT_BH(config,subDir);        
%         
%         fileID      = fopen(filename,'a');
%         fprintf(fileID,'# %s \n','** 2nd Iteration **');
%         fprintf(fileID,'# %s \n',['maxComp_y2: ',num2str(config.optsNum.maxComp_y2,3)]);
%         fprintf(fileID,'# %s \n',['CA Grid: ',num2str(config.optsNum.PhysArea.alpha_deg,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
%         fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
%         fclose(fileID);
%         
%          %(3) 3rd Iteration        
%          config.optsNum.PhysArea.alpha_deg = round(sol.Measured_StaticCA*180/pi);
         config.optsNum.maxComp_y2         = 25; 
         
         sol = ThreePhaseContactLine_FMT_BH(config,subDir);
         
%         config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;      
%         degAngle                          = sol.Measured_StaticCA*180/pi;
%         
         fileID      = fopen(filename,'a');
         fprintf(fileID,'# %s \n','** 3rd Iteration **');
         fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
         fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
         fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
         fclose(fileID);
        
    end
end



%         PhysArea = struct('N',[40,40],'L1_Skewed',2,'L2',2,'y2wall',0.,...
%                           'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);
% 
%         PhysArea.Conv  = struct('L',1,'L2',[],'N',[20,20]);
% 
%         Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
%                            'Ncircle',1,'N1disc',35,'N2disc',34);
% 
%         optsNum = struct('PhysArea',PhysArea,...
%                          'FexNum',Fex_Num,...
%                          'maxComp_y2',10,...
%                          'y1Shift',0,...c 
%                          'plotTimes_T',100,...
%                          'plotTimes_Interval',0.1);
% 
%         V1 = struct('V1DV1','Vext_Cart_7',...
%                             'epsilon_w',0.74,... %0.482218,...
%                             'epsilon_w_max',0.3,....
%                             'tau',20);
%                         
%         V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'R',1);
% 
%         optsPhys = struct('V1',V1,'V2',V2,...                      
%                           'kBT',0.75,...   
%                           'gamma',2,...
%                           'Inertial',1,...
%                           'Dmu',0.0,'nSpecies',1,...
%                           'sigmaS',1);      
%                       
%         configuration = v2struct(optsNum,optsPhys);
%         configName    = SaveConfig(configuration,'Configurations');