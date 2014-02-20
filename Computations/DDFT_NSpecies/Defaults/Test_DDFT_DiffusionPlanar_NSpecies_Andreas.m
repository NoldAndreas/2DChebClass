function [optsNum,optsPhys] = Test_DDFT_DiffusionPlanar_NSpecies()

    %Numerical Parameters    
    Phys_Area = struct('N1',20,'N2',20, ...
                       'L1',2,'L2',2);
    
    Plot_Area = struct('y1Min',-10,'y1Max',10,'N1',100,...
                       'y2Min',-10,'y2Max',10,'N2',100);
    
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
                   
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_DiffusionPlanar_NSpecies',...
                     'plotTimes',0:7/100:7, ...
                     'FexNum',FexNum);    
                                  
    epsilonS=0*[ 1 1 1 ;
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
    
    V1       = struct('V1DV1','quadBump',...
                      'V0',2,'y10',2,'y20',2,'tau',10);                                         
                  
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);        
        
    optsPhys = struct('V1',V1,'V2',V2,...
                      'kBT',0.7,...
                      'nParticlesS',[20;20;20]); 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum,optsPlot);                 

end                 

