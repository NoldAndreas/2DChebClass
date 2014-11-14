function EX = DDFT_DiffusionInfCapillary_FMT_Case1(doHI)

    if(nargin==0)
        doHI = false;
    end

    Phys_Area = struct('shape','InfCapillary_FMT','N',[60;40],... %
                       'L1',3,'y2Min',0,'y2Max',4,'N2bound',10);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',80,...
                       'y2Min',Phys_Area.y2Min,...
                       'y2Max',Phys_Area.y2Max,'N2',80);

    Sub_Area = struct('shape','Box','y1Min',-1,'y1Max',1,'N',[40,40],...
                      'y2Min',0,'y2Max',4);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',34,'N2disc',34);
	V2Num   = struct('Fex','SplitDisk','L',1,'L1',1,'L2',[],'N',[34,34]);            

    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'SubArea',Sub_Area,...
                     'FexNum',Fex_Num,...         
                     'plotTimes',0:0.5:6,...
                     'V2Num',V2Num);
    
    %V1       = struct('V1DV1','Vext_Cart_1','V0',0.2,...
    %                    'grav',1,'y10',0,'y20',6);   

% V1 = struct('V1DV1','Vext_Cart_Slit_Static','epsilon_w',[1 1 1 1]);
    %V1 = struct('V1DV1','Vext_BarkerHenderson_HardWall','epsilon_w',1.0);
   % V1 = struct('V1DV1','zeroPotential');      

    V1       = struct('V1DV1','Vext_Cart_2','V0',0.2,...
                        'grav',1,'y10',0,'y20',6);  
    
    V2 = struct('V2DV2','BarkerHenderson_2D','epsilon',1,'LJsigma',1);                     
    
    optsPhys = struct('V1',V1,'V2',V2,...    
                      'kBT',0.75,...
                      'nParticlesS',5,...
                      'nSpecies',1,'sigmaS',1); 

    if(doHI)
        optsPhys.HI = HI;
        optsNum.HINum = HI_Num;
    end
                  
    optsPlot.doDDFTPlots=true;
    
    AddPaths();
    EX   = DDFT_2D(v2struct(optsPhys,optsNum));
    EX.Preprocess();
    EX.ComputeEquilibrium();%EX.optsPhys.rhoGas_sat);
    
    EX.IDC.plot(EX.GetRhoEq,'SC');
    EX.ComputeDynamics();
    EX.PlotDynamics();
end                 
