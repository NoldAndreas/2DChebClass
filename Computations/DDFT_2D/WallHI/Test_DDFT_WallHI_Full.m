function EX = Test_DDFT_WallHI_Full(HItype)

    if(nargin==0 || strcmp(HItype,'none'))
        doHI = false;
        HItype = '';
    else
        doHI = true;
    end

    Phys_Area = struct('shape','HalfSpace_FMT','N',[25;25],'L1',2,'L2',2, ...
                       'y2wall',0,'N2bound',10,'h',1,'L2_AD',1,'alpha_deg',90);                    
    
    Sub_Area = struct('shape','Box','y1Min',-5,'y1Max',5,'N',[20,20],...
                      'y2Min',0.5,'y2Max',1);                   
                   
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',5,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);

    if(doHI)
        HItype = [HItype '_2D_noConv'];
    end
                   
    %HItype = 'Oseen_2D_noConv';
    %HItype = 'TwoOseen_2D_noConv';
    %HItype = 'FullWallHI_2D_noConv';
                                
    HI_Num    = struct('N',[30;30],'L',2,'HI11','noHI_2D','HI12',HItype, ...
                      'HIPreprocess', 'RotnePragerPreprocess2D',...
                      'HIWallFull',true,'doConv',false,...
                      'Wall','SelfWallTerm');


    tMax = 0.3;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'SubArea',Sub_Area,...
                     'FexNum',Fex_Num,...                    
                     'plotTimes',0:tMax/100:tMax);
    
	sigmaS  = 1;
    sigmaHS = 1;

%     V1       = struct('V1DV1','V1_Well_gravity',...
%                       'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
%                       'y10',0,'y20',2,'g',1,'gcut',10); 
    
    V1       = struct('V1DV1','V1_Well_gravity',...
                      'V0',0.01,'V0add',3,'tau',0.1,'sigma1Add',0.5,'sigma2Add',0.5, ...
                      'y10',0,'y20',1,'g',-1,'gcut',10); 




    HI       = struct('sigmaS',sigmaS,'sigmaHS',sigmaHS);
    
    optsPhys = struct('V1',V1,  ...                                            
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',20,'sigmaS',sigmaS);

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;

    EX = DDFTDynamics(optsPhys,optsNum,optsPlot);
    
end                 


