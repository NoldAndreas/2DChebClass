function [optsNum,optsPhys] = Test_DDFT_EqBox_2Phase_2()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-10,'y1Max',10,'L1',5,'N1',20,...
                       'y2Min',-4,'y2Max',4,'N2',20);
    
    Plot_Area = Phys_Area;  Plot_Area.N1 = 100; Plot_Area.N2 = 100;
   
            
    optsNum = struct('PhysArea',Phys_Area,'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_EqBox_2Phase_2','Cutoff',5);                            
                      
    optsPhys = struct('V1DV1','Vext_Cart_Capillary_Static',...                      
                      'V0',0.0,'y10',-5,'y20',0,...
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,...
                      'V2DV2','Phi2DLongRange');                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

