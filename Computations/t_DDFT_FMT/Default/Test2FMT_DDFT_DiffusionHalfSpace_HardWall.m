function [optsNum,optsPhys,optsPlot] = Test2FMT_DDFT_DiffusionHalfSpace_HardWall()

    %This is equivalent to:

    %Numerical Parameters    
    Phys_Area = struct('N',[60;60],'L1',2,'L2',1.,'y2wall',0.,...
                       'N2bound',16,'h',1,'L2_AD',1.);
    
    Plot_Area = struct('y1Min',-2,'y1Max',2,'N1',100,...
                       'y2Min',0.5,'y2Max',3,'N2',100);
                   
    Fex_Num   = struct('Fex','FMTRosenfeld_3DFluid',...
                       'Ncircle',1,'N1disc',24,'N2disc',24);
    
    %Sub_Area  = Phys_Area;
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,...
                     'FexNum',Fex_Num,...
                     'DDFTCode','CheckFMTComputation',...%DDFT_DiffusionHalfSpace_NSpeciesHardWall',...
                     'Tmax',7,'TN',50,...
                     'name','default',...
                     'Accuracy_Averaging',1e-6);
                     
%    FexMatrices = {'Polar_SpectralFourier_FMTMatricesFull', ...
%               'Polar_SpectralFourier_FMTMatricesFull_Roth', ...
%               'Polar_SpectralFourier_FMTMatricesFull_3D'};
%    Fex         = {'Polar_SpectralFourier_FMT', ...
%               'Polar_SpectralFourier_FMT_Roth', ...
%               'Polar_SpectralFourier_FMT_3D'};                                  

    sigmaS = [1];%   1.1 ;
              %1.1 1.2 ];        
              
    V1       = struct('V1DV1','Vext_Cart_6','V0',0,...
                      'grav_start',-2,'grav_end',0,'T',5,...
                      'L1',1  ,'L2',1);              
        
    %V1       = struct('V1DV1','zeroPotential',...  %Vext_Cart_6
    %                  'V0',0.1,'V0r',1,'alphar',1,'tau',1,'rV',1);                                         
                  
    V2       = struct('V2DV2','zeroPotential');                                 
    
    optsPhys = struct('V1',V1,'V2',V2,...                                            
                      'kBT',1,'eta',0.4783,...%4783,...%0.4257,...
                      'nParticlesS',10,...
                      'sigmaS',sigmaS); %[10;10]
                 
    lineColourDDFT={{'r','b','g'}};            
    optsPlot = struct('lineColourDDFT',lineColourDDFT);
    optsPlot.doDDFTPlots=true;
                  
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    
    if(nargout == 0)
        f(optsPhys,optsNum,optsPlot);                 
    end
%     for l2 = 5.0:0.5:6
%         optsNum.PhysArea.L2 = l2;
%         for n = 80:1:80
%             optsNum.PhysArea.N(2) = n;
%             name = ['HardWall_N=',num2str(n),'L2=',num2str(l2)];
%             f(optsPhys,optsNum,optsPlot,name);                 
%         end
%     end

end                 

