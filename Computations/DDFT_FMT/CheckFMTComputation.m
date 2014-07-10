  function CheckFMTComputation()
%************************************************************************* 
% define
%   mu_s_i := kBT*log(rho_i) + sum( int(rho_j(r')*Phi2D(r-r'),dr') , 
%               j = 1..N) +mu_HS(rho_i) + V_ext_i - mu_i
% Equilibrium, (solve for each species i)
%  (EQ i,1) mu_s_i     = 0
%  (EQ i,2) int(rho_i) = NoParticles_i
% Dynamics: 
%*************************************************************************   

    %********************************************
    %***************  Parameters ****************
    %********************************************
    disp('** CheckFMTComputationSkewed **');
    %Numerical Parameters    
    Phys_Area = struct('N',[1;150],'L1',2,'L2',2.,'y2wall',0.,...
                       'N2bound',30,'h',1,'L2_AD',2.,...
                       'alpha',pi/2);

    Plot_Area = struct('y1Min',-5,'y1Max',5,'N1',100,...
                       'y2Min',0.5,'y2Max',6,'N2',100);

    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',30,'N2disc',30);

    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num);

    V1       = struct('V1DV1','zeroPotential');              
    V2       = struct('V2DV2','zeroPotential');                                 

    optsPhys = struct('V1',V1,'V2',V2,...                                            
                      'kBT',1,'eta',0.4257,...
                      'sigmaS',1);%0.4783,...%0.4257,...                      

    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;

    AddPaths();    
    close all;  

    %************************************************
    %***************  Initialization ****************
    %************************************************
    PhysArea  = optsNum.PhysArea;           
    Rs        = diag(optsPhys.sigmaS)/2;     
    HS        = HalfSpace_FMT(PhysArea,Rs);
    
    %************************************************
    %****************  Preprocess  ******************
    %************************************************
    
    tic
    fprintf(1,'Computing Fex matrices ...\n');   
    
    params         = optsPhys;
    params         = rmfield(params,'V1');
    params.FexNum  = optsNum.FexNum;
    params.PhysArea = optsNum.PhysArea;
    
    params.Pts     = HS.Pts;     
    params.Polar   = 'cart';      
    func           = str2func(['FexMatrices_',optsNum.FexNum.Fex]);            
    params         = rmfield(params,'eta');
    
    IntMatrFex_2D  = DataStorage(['HalfSpace_FMT' filesep func2str(func)],func,params,HS); %true 
    
    CheckAverageDensities_Rosenfeld_3D(HS,IntMatrFex_2D);
    
    optsPhys.rho_iguess = optsPhys.eta*6/pi;
    rho_ic1D = FMT_1D_HardWall(HS,IntMatrFex_2D,optsPhys,optsNum);        
    
    

end