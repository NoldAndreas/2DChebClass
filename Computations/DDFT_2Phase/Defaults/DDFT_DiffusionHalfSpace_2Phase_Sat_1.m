function [optsNum,optsPhys] = DDFT_DiffusionHalfSpace_2Phase_Sat_1()

    %Numerical Parameters    
    Phys_Area = struct('N',[10,50],'y2Min',0,'L1',3,'L2',5);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,'N2',100,...
                       'y2Min',0,'y2Max',10);
    
    %Sub_Area  = Phys_Area;
    Sub_Area = struct('y1Min',-1,'y1Max',1,'N',[20,20],...
                      'y2Min',0,'y2Max',2);    
                  
	Conv      = struct('L1',2,'L2',1.,'N',[40,40]);
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'plotTimes',0:0.1:6,...
                     'Conv',Conv,...
                     'DDFTCode','DDFT_DiffusionHalfSpace_2Phase_Sat');
                      

    V1 = struct('V1DV1','Vext_Cart_5',...
                      'V0',0.0,'epsilon_w',1,'y10',2,'y20',0,'tau',1,...
                      'epsilon_w_end',0.0,'wallAxis',2);
                  
    V2 = struct('V2DV2','Phi2DLongRange','epsilon',1);
                 
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','CarnahanStarling',...
                      'kBT',0.7,...                      
                      'Dmu',-0.05,'nSpecies',1);

    %Run File
    if(nargout == 0)
        f = str2func(optsNum.DDFTCode);    
        f(optsPhys,optsNum);
    end
end                 

