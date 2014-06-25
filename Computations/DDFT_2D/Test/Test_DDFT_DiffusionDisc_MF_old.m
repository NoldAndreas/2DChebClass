function [output,optsNum,optsPhys,optsPlot] = Test_DDFT_DiffusionDisc_MF_old()

    Phys_Area = struct('shape','Disc','N',[20;20],'R',4);
    
    Plot_Area = struct('y1Min',0,'y1Max',2,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
                   
    Sub_Area = struct('shape','Box','N',[20,20],...
                      'y1Min',-1,'y1Max',1,...
                      'y2Min',-1,'y2Max',1);  
                   
    Fex_Num    = struct('Fex','Meanfield','N',[20,20],'L',1);
%    Fex_Num   = struct('Fex','FMTRosenfeld',...
%                       'Ncircle',10,'N1disc',10,'N2disc',10);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'SubArea',Sub_Area,...
                     'FexNum',Fex_Num,'DDFTCode','DDFT_Diffusion_2D',...                     
                     'plotTimes',0:7/100:7);

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
        
        
    V1 = struct('V1DV1','V1_Test_Disc','V0',0.5,'grav',1); 
                           
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);
    
    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',[20;20;20]); 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    
    config = v2struct(optsPhys,optsNum);
    
    %****************************
    AddPaths();
%     EX     = DDFT_2D(config);
%     EX.Preprocess();
%     EX.ComputeEquilibrium();
    
    
    f = str2func(optsNum.DDFTCode);
    output = f(optsPhys,optsNum);                 

end                 

