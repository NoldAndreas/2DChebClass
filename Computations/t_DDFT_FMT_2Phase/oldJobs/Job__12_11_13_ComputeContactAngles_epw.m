function Job__12_11_13_ComputeContactAngles_epw()
%This Job should compute the contact lines for a range of contact angles, 
% such that a full study of e.g. the line tensions is made possible
% For enhanced accuracy, 60x60 points are used. 
% No dynamics is to be computed
    global dirData
    AddPaths();

    PhysArea = struct('N',[50,50],'L1',2,'L2',2,'y2wall',0.,...
                      'N2bound',24,'h',1,'L2_AD',2.,'alpha_deg',90);

    PhysArea.Conv  = struct('L',1,'L2',[],'N',[20,20]);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',35,'N2disc',34);

    optsNum = struct('PhysArea',PhysArea,...
                     'FexNum',Fex_Num,...
                     'maxComp_y2',[],...
                     'y1Shift',0);

	%For Vext_BarkerHenderson_HardWall, at kbT = 0.75, => 
    % epsilon_w  | CA
    %----------------------
    %  1.0       | 90 deg
    %  1.2       | 67 deg
    %  1.4       | 38.6 deg
    %  1.5       | 15 deg
    %-----------------------
    V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.0);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'R',1); 

    optsPhys = struct('V1',V1,'V2',V2,...                   
                      'kBT',0.75,...                                                    
                      'Dmu',0.0,'nSpecies',1,...
                      'sigmaS',1);      

    config = v2struct(optsNum,optsPhys);                    
    config.optsNum.PhysArea.alpha_deg = 90;
    
    %***********************************************************
    %Setup result file for this Job
    subDir      = 'FMT_ContactLine_Equilibrium_BH_40X40_epw';
    
    
    time = clock();
    filename   = ([dirData filesep subDir filesep]);
    filename   = [filename,'Job__12_11_13_ComputeContactAngles_epw_',...
               num2str(time(1)),'_'... % Returns year as character
               num2str(time(2)),'_'... % Returns month as character
               num2str(time(3)),'_'... % Returns day as char
               num2str(time(4)),'_'... % returns hour as char..
               num2str(time(5)),...    %returns minute as char                    
               '.txt'];
           
  %  filename    = [dirData filesep subDir filesep 'Job__12_11_13_ComputeContactAngles_epw.txt'];
    Struct2File(filename,config,['Computed at ',datestr(now)]);
    
    epw    = [1.0:0.1:1.4];%,1.11,1.12];
    
    degAngle    = 90;
    %epw       = [1.11,1.12];
    %degAngle  = 22.0149;
    
    for i = 1:length(epw)
        
        config.optsPhys.V1.epsilon_w  = epw(i);
        fileID = fopen(filename,'a');
        fprintf(fileID,'# ******************* \n');
        fprintf(fileID,'# %s \n',['epw =  ',num2str(epw(i),3)]);
        fclose(fileID);
        
        config.optsNum.maxComp_y2         = -1; %No 2D Computation
        config.optsNum.PhysArea.alpha_deg = 90; 
        sol = ThreePhaseContactLine_FMT_BH(config,subDir);
        fileID      = fopen(filename,'a');
        fprintf(fileID,'# %s \n',['Young CA for 90 deg grid: ',num2str(sol.Young_CA*180/pi,3)]);
        fclose(fileID);
        
        config.optsNum.PhysArea.alpha_deg = degAngle;         
        
        config.optsNum.maxComp_y2         = 10;
        sol = ThreePhaseContactLine_FMT_BH(config,subDir);                
        config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi; 
        
        fileID      = fopen(filename,'a');
        fprintf(fileID,'# %s \n','** 1st Iteration **');
        fprintf(fileID,'# %s \n',['CA Grid: ',num2str(config.optsNum.PhysArea.alpha_deg,3),'[deg]']);
        fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
        fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
        fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
        fclose(fileID);
        
%         config.optsNum.maxComp_y2         = 15;
%         sol = ThreePhaseContactLine_FMT_BH(config,subDir);
%         config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;                       
%         
%         fileID      = fopen(filename,'a');
%         fprintf(fileID,'# %s \n','** 2nd Iteration **');
%         fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
%         fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
%         fclose(fileID);
%         
%         config.optsNum.maxComp_y2         = 20;
%         sol = ThreePhaseContactLine_FMT_BH(config,subDir);
%         config.optsNum.PhysArea.alpha_deg = sol.Measured_StaticCA*180/pi;      
%         degAngle = sol.Measured_StaticCA*180/pi;
%         
%         fileID      = fopen(filename,'a');
%         fprintf(fileID,'# %s \n','** 3rd Iteration **');
%         fprintf(fileID,'# %s \n',['CA from Young Eq: ',num2str(sol.Young_CA*180/pi,3),'[deg]']);
%         fprintf(fileID,'# %s \n',['CA measured: ',num2str(sol.Measured_StaticCA*180/pi,3),' [deg]']);
%         fprintf(fileID,'# %s \n',['Data saved in ',sol.Filename]);
%         fclose(fileID);
        
    end
end