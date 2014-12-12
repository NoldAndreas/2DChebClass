function [EX,res] = Test_DDFT_InertiaInfDisc_MF_osc(inertial)

    if(nargin==0)
        inertial = true;
    end

    Phys_Area = struct('shape','InfDisc','N',[20;10], ...
                       'y1Min',0,'y1Max',inf,'L',4,...
                       'y2Min',0,'y2Max',2*pi);
    
    Plot_Area = struct('y1Min',0,'y1Max',5,'N1',100,...
                       'y2Min',0,'y2Max',2*pi,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    gamma = 2;
    tMax = pi*gamma/2;
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',FexNum,...                     
                     'plotTimes',0:tMax/100:tMax); 

    epsilonS=2;
    alphaS=1;
            
    V1       = struct('V1DV1','V1_Well_osc',...
                      'V0',0.0001,'V0add',2,'sigma1Add',5,'sigma2Add',5,'y10',0,'y20',0,'tau',tMax/pi); 
                           
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);
    
    optsPhys = struct('Inertial',inertial,'V1',V1,'V2',V2,...
                      'kBT',1,'mS',1,'gammaS',gamma, ...
                      'nParticlesS',10); 
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
      
    [EX,res] = DDFTDynamics(optsPhys,optsNum,optsPlot);
end                 

