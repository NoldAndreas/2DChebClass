function [EX,res] = DDFT_DiffusionInfSpace_FMT(doHI)

    if(nargin==0)
        doHI = true;
    end

    Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[30,30],'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld',...
                       'Ncircle',10,'N1disc',10,'N2disc',10);
                   
    HI_Num    = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
                      'HIPreprocess', 'RotnePragerPreprocess2D');  
    
    tMax = 0.25;    
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...                            
                     'plotTimes',0:tMax/100:tMax);%'DDFTCode','DDFT_Diffusion_2D',...

    sigmaS  = 1;
    sigmaHS = 0.5;

    V1       = struct('V1DV1','V1_Triangle',...
                      'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
                       'y10',-1,'y20',-0.5,'y11',1,'y21',0,'y12',0,'y22',0.5);                      
    
    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1,...                                            
                      'kBT',1,'mS',1,'gammaS',1,'sigmaS',sigmaS,...
                      'nParticlesS',20); 

    optsPlot.lineColourDDFT  = {'r'};
    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
        optsPlot.lineColourDDFT  = {'b'};
    end
                  
    optsPlot.doDDFTPlots=true;
    
    [EX,res] = DDFTDynamics(optsPhys,optsNum,optsPlot);

end                 

