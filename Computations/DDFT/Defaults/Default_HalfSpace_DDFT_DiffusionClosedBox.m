function [optsNum,optsPhys] = Default_HalfSpace_DDFT_DiffusionClosedBox()
    %Numerical Parameters    
    Phys_Area = struct('y1Min',-inf,'y1Max',inf,'L1',3,'L2',2,'N',[20;20],...
                       'y2Min',0);
                   
    Plot_Area       = struct('y1Min',-5,'y1Max',5,'L1',3,'L2',2,...
                       'N2',100,'N1',100,'y2Min',0,'y2Max',4);    
    
    Sub_Area        = struct('y1Min',-3,'y1Max',3,'N',[20;20],...
                             'y2Min',-2,'y2Max',2);
    
    Conv      = struct('L',2,'L2',2,'N',[30,30]);%'ep2Conv',0.1
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'TN',50,'Tmax',5,...
                     'Conv',Conv,...
                     'DDFTCode','DDFT_DiffusionInfCapillary');                            
                 
    V1 = struct('V1DV1','Vext_Cart_6','V0',0.,'grav_start',-1.0,'grav_end',0,'T',2);
    V2 = struct('V2DV2','Gaussian','alpha',2,'epsilon',0.06);
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                     'kBT',0.7,'mu',0,'nSpecies',1);

    if(nargout == 0)
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        f(optsPhys,optsNum);    
    end

end       