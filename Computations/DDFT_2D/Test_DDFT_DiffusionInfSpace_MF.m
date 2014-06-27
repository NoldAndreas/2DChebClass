function [output,optsNum,optsPhys,optsPlot] = Test_DDFT_DiffusionInfSpace_MF()

    Phys_Area = struct('shape','InfSpace','N',[20;20], ...
                       'y1Min',-inf,'y1Max',inf,'L1',4,...
                       'y2Min',-inf,'y2Max',inf,'L2',4);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,...
                       'y2Min',-10,'y2Max',10,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',FexNum,...                     
                     'plotTimes',0:4/100:4); 

    epsilonS=2*[ 1 1 1 ;
    1 1 1 ;
    1 1 1 ] ;

    alpha1 = 0.5;
    alpha2 = 1;
    alpha3 = 1.5;

    alpha12=(alpha1+alpha2)/2;
    alpha13=(alpha1+alpha3)/2;
    alpha23=(alpha2+alpha3)/2;

    alphaS=[alpha1  alpha12 alpha13 ; 
            alpha12 alpha2  alpha23 ;
            alpha13 alpha23 alpha3  ];
        
        
    V1       = struct('V1DV1','V1_Rotating_Cart',...
                      'V0',0.05,'V0r',1,'alphar',1,'tau',1,'rV',1); 
                           
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);
    
    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',[10;10;10]); 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
      
	config = v2struct(optsPhys,optsNum);
    
    AddPaths();
    EX     = DDFT_2D(config);
    EX.Preprocess();
    EX.ComputeEquilibrium();
    EX.ComputeDynamics();    
end                 
