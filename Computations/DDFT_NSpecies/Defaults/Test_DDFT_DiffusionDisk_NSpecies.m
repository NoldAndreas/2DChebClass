 function [optsNum,optsPhys] = Test_DDFT_DiffusionDisk_NSpecies()

    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',10,'N',[20,20],...
                       'y2Min',0,'y2Max',2*pi);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_DiffusionPolar_NSpecies',...
                     'plotTimes',0:7/100:7,...
                     'FexNum',FexNum);    
                                  
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
    
    potParams2Names ={{'epsilon','alpha'}};        
    
    optsPhys = struct('V1DV1','quadSwitchDisk',...
                      'V2DV2','Gaussian',...
                      'potParams2Names',potParams2Names, ...
                      'epsilonS',epsilonS,'alphaS',alphaS, ...
                      'V0',2,'y10',0,'y20',0,'tau',10,...
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

