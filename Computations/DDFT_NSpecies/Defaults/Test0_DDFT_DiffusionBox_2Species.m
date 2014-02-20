function [optsNum,optsPhys] = Test0_DDFT_DiffusionBox_2Species(exec)

    N1 = 20; N2 = 20;
    %Numerical Parameters    
    Phys_Area = struct('y1Min',0,'y1Max',10,'N',[N1,N2],...
                       'y2Min',0,'y2Max',10);
    
    Plot_Area = struct('y1Min',0,'y1Max',10,'N1',100,...
                       'y2Min',0,'y2Max',10,'N2',100);
                   
    FexNum = struct('Fex','Meanfield');   
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'DDFTCode','DDFT_DiffusionBox_NSpecies',...
                     'plotTimes',0:1:70,...
                     'FexNum',FexNum);    
               
    epsilonS=  [ 1   -0.1 ;
                 -0.1  1 ] ;

    epsilon_w1     = repmat([1,1],N1*N2,1);
    epsilon_w1_end = repmat([1,0],N1*N2,1);
    epsilon_w2     = repmat([0,0],N1*N2,1);
    epsilon_w2_end = repmat([0,1],N1*N2,1);                
    
%    V1       = struct('V1DV1','rotating2','V0',0.05,'V0r',1,'alphar',1,'tau',1,'rV',1);                                                           
%    V2       = struct('V2DV2','hardSphere','sigmaS',sigmaS);                                 

    V1        = struct('V1DV1','Vext_Cart_2Species_1',...
                      'V0',0.0,'y10',4,'y20',4,'tau',5,...
                      'epsilon_w1',epsilon_w1,'epsilon_w1_end',epsilon_w1_end,...
                      'epsilon_w2',epsilon_w2,'epsilon_w2_end',epsilon_w2_end);    
    V2       = struct('V2DV2','Phi2DLongRange','epsilon',epsilonS);                                     
                      
    optsPhys = struct('V1',V1,'V2',V2,...
                      'HSBulk','MuCarnahanStarling',...
                     'kBT',0.7,...
                     'nParticlesS',[30;30]);
                 
    %Run File
    if((nargin == 0) || ((nargin == 1) && exec))
        AddPaths();
        f = str2func(optsNum.DDFTCode);
        f(optsPhys,optsNum);                 
    end

end                 

