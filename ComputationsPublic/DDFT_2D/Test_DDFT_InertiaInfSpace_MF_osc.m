function [EX,res] = Test_DDFT_InertiaInfSpace_MF_osc(inertial)

    if(nargin==0)
        inertial = true;
    end

    Phys_Area = struct('shape','InfSpace','N',[30;30], ...
                       'y1Min',-inf,'y1Max',inf,'L1',3,...
                       'y2Min',-inf,'y2Max',inf,'L2',3);
    
    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',-5,'y2Max',5,'N2',100);
                   
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

