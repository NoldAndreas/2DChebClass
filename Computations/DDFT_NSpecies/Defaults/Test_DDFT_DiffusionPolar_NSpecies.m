function [optsNum,optsPhys,optsPlot] = Test_DDFT_DiffusionPolar_NSpecies()

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',inf,'N',[15,10],'L1',4,...
                       'y2Min',0,'y2Max',2*pi);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',FexNum,...
                     'DDFTCode','DDFT_DiffusionPolar_NSpecies',...
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
        
    potParams2Names = {{'epsilon','alpha'}};
        
    V1       = struct('V1DV1','rotating2',...
                      'V0',0.05,'V0r',1,'alphar',1,'tau',1,'rV',1);                                         
                  
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);
    
    optsPhys = struct('V1',V1,'V2',V2,...                         
                      'kBT',1,...
                      'potParams2Names',potParams2Names,...
                      'nParticlesS',[10;10;10]); 
                  
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum,optsPlot);                 

end                 

