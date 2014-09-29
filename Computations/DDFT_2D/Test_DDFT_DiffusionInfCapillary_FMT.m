function EX = Test_DDFT_DiffusionInfCapillary_FMT(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','InfCapillary_FMT','N',[1;20],... %
                       'L1',2,'y2Min',0,'y2Max',8,'N2bound',10);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',80,...
                       'y2Min',Phys_Area.y2Min,...
                       'y2Max',Phys_Area.y2Max,'N2',80);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',20,'N1disc',20,'N2disc',20);
	V2Num   = struct('Fex','SplitDisk','L',1,'L1',1,'L2',[],'N',[10,10]);    
    
    
    tMax = 10;

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...         
                     'V2Num',V2Num,...
                     'plotTimes',0:tMax/100:tMax);                     
    
%    V1       = struct('V1DV1','Vext_Cart_1','V0',0.2,...
%                        'grav',1,'y10',0,'y20',6);                        	                        
    V1 = struct('V1DV1','Vext_Cart_Slit_Static','epsilon_w',[1 1 2 2]);
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1);                     
    
    optsPhys = struct('V1',V1,'V2',V2,...    
                      'kBT',0.75,'Dmu',0.0,'nSpecies',1,'sigmaS',1); 

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;
    
    EX = DDFTDynamics(optsPhys,optsNum,optsPlot);
    
end                 
