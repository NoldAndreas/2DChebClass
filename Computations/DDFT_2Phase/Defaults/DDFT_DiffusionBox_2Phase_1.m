                 
function [optsNum,optsPhys] = DDFT_DiffusionBox_2Phase_1()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',10,'N',[25,25],...
                       'y2Min',-6,'y2Max',6);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',-6,'y2Max',6,'L2',10,'N2',100);
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',5,'y1Max',10,'N',[20,20],...
                      'y2Min',0,'y2Max',5);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'DDFTCode','DDFT_DiffusionBox_2Phase',...
                     'plotTimes',0:1:20);                            
                     
    V1 = struct('V1DV1','Vext_Cart_5',...
                      'V0',0.02,'epsilon_w',1,'epsilon_w_end',0,...
                      'y10',2,'y20',0,'tau',1);
                  
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,'nParticlesS',40);
                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

