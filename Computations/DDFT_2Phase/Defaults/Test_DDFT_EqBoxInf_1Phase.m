function [optsNum,optsPhys] = Test_DDFT_EqBoxInf_1Phase()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'L1',5,'N1',25,...
                       'y2Min',-4,'y2Max',4,'N2',20);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'L1',5,'N1',100,...
                       'y2Min',-4,'y2Max',4,'N2',100);
            
    optsNum = struct('PhysArea',Phys_Area,'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_EqBoxInf_1Phase','Cutoff',5);                            
                      
    optsPhys = struct('V1DV1','Vext_Cart_Slit_Static',...                                            
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,...
                      'V2DV2','Phi2DLongRange',...
                      'V2DV2_1D','Phi1DLongRange');                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

