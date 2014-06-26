function Test1_DFT_1D_LiqGas_2Phase()

    %Numerical Parameters    
    Phys_Area = struct('yMin',-inf,'yMax',inf,'L',5,'N',21);    
    Plot_Area = struct('yMin',-20,'yMax',20,'N',500);    
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'DDFTCode','DFT_1D_LiqGas_2Phase',...
                     'plotTimes',0:5:500,...
                     'Cutoff',4,'Conv_N',40);                            
                      
    optsPhys = struct('HSBulk','CarnahanStarling','kBT',0.7,...
                     'V2DV2','Phi1DLongRange');
                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);

end                 

