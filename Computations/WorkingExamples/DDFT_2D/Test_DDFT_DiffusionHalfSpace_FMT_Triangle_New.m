function [EX,res] = Test_DDFT_DiffusionHalfSpace_FMT_Triangle_New(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','HalfSpace_FMT','N',[20;20],'L1',2,'L2',2, ...
                       'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90);                    
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

    %HItype = 'FullWallHI_2D_noConv';
    %HItype = 'Oseen_2D_noConv';
    HItype = 'TwoOseen_2D_noConv';

    HI_Num    = struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12',HItype, ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false);
    
%     HI_Num    = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D_noConv', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D',...
%                       'HIWallFull',true,'doConv',false);  
 
%     HI_Num    = struct('N',[20;20],'L',2,'HI11','noHI_2D','HI12','RP12_2D', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D',...
%                       'HIWallFull',true,'doConv',true); 

%     HI_Num    = struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12','Oseen_2D_noConv', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D',...
%                       'HIWallFull',true,'doConv',false);

%     HI_Num    = struct('N',[10;10],'L',2,'HI11','noHI_2D','HI12','TwoOseen_2D_noConv', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D',...
%                       'HIWallFull',true,'doConv',false);

%     HI_Num    = struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12','FullWallHI_2D_noConv', ...
%                       'HIPreprocess', 'RotnePragerPreprocess2D',...
%                       'HIWallFull',true,'doConv',false);


    tMax = 0.6;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...                    
                     'plotTimes',0:tMax/100:tMax);
    
	sigmaS  = 1;
    %sigmaHS = 0.5;
    sigmaHS = 1;
    
    V1       = struct('V1DV1','V1_Triangle',...
                      'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
                      'y10',-1,'y20',2.5,'y11',1,'y21',3,'y12',0,'y22',3.5); 

    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1,  ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',20,'sigmaS',sigmaS);

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;

    [EX,res] = DDFTDynamics(optsPhys,optsNum,optsPlot);
    
end                 


