function [optsNum,optsPhys] = DDFT_DiffusionBox_2Phase_Sat_1()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',10,'N',[20,30],...
                       'y2Min',-10,'y2Max',10);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,'N2',100,...
                       'y2Min',-10,'y2Max',10);
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',0,'y1Max',5,'N',[20,20],...
                      'y2Min',-10,'y2Max',10);    
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.2:6,...
                     'DDFTCode','DDFT_DiffusionBox_2Phase_Sat');                            
                      

    V1 = struct('V1DV1','Vext_Cart_5',...
                      'V0',0.0,'epsilon_w',1,'y10',2,'y20',0,'tau',1,...
                      'epsilon_w_end',0.0);
                  
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,...                      
                      'nParticlesS',40);

    %Run File
    if(nargout == 0)
        f = str2func(optsNum.DDFTCode);    
        f(optsPhys,optsNum);
    end
end                 

