function [optsNum,optsPhys] = Test0_DDFT_DiffusionDisk_2Species()

    N1 = 20; N2 = 20;
    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',4,'L1',2,'N',[N1,N2],...
                       'y2Min',0,'y2Max',2*pi,'L2',5);
    
    Plot_Area = struct('y1Min',0,'y1Max',4,'L1',2,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'L2',5,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[10,10],'L',1);                   
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_DiffusionPolar_NSpecies',...
                     'plotTimes',0:0.2:5,...
                     'FexNum',FexNum); 
                 
    epsilonS         =  [ 1    -0.1;
                          -0.1  1];
    grav             = repmat([1,-1],N1*N2,1);
                     
    V1               = struct('V1DV1','Vext_Pol_2Species_1',...
                            'V0',0.3,'tau',5,'grav',grav);
    V2               = struct('V2DV2','Phi2DLongRange','epsilon',epsilonS);
                      
    optsPhys = struct('V1',V1,'V2',V2,...                                            
                      'HSBulk','MuCarnahanStarling',...
                      'kBT',0.7,...
                      'nParticlesS',[10;10]);
                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

