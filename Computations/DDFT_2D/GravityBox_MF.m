 function [EX,res] = GravityBox_MF(theta)

    if(nargin==0)
        theta = 0;
    end
 
    Phys_Area = struct('shape','Box','N',[20;20], ...
                       'y1Min',0,'y1Max',10,...
                       'y2Min',0,'y2Max',20);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',0,'y2Max',20,'N2',100);
                   
    FexNum = struct('Fex','Meanfield','N',[20,20],'L',1);
    
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'V2Num',FexNum,...                     
                     'plotTimes',0:7/100:7);

    epsilonS = 0;

    alphaS = 1;    
        
    V1       = struct('V1DV1','gravity','g',1,'theta',theta); 
                           
    V2       = struct('V2DV2','Gaussian','epsilon',epsilonS,'alpha',alphaS);
    
    optsPhys = struct('Inertial',false,'V1',V1,'V2',V2,...
                      'kBT',1,'mS',1,'gammaS',1, ...
                      'nParticlesS',20); 

    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    [EX,res] = DDFTDynamics(optsPhys,optsNum,optsPlot);
end                 

