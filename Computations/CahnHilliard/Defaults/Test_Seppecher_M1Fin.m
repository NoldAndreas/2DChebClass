                 
function [optsNum,optsPhys] = Test_Seppecher_M1Fin()

    %Numerical Parameters    
    Phys_Area = struct('y2Min',0,'y2Max',10,'N2',40,'L1',2,...
                       'y1Min',-10,'y1Max',10,'N1',30);
    
    Plot_Area = Phys_Area;
    Plot_Area.N1 = 100;  Plot_Area.N2 = 150;    
                        
    Sub_Area = struct('y2Min',0,'y2Max',5,'N2',30,'L1',2,...
                       'y1Min',0,'y1Max',5,'N1',30);
        
    optsNum   = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'DDFTCode','Seppecher_M1Fin',...
                     'plotTimes',0:0.1:5);                            
                      
    optsPhys  = struct('theta',pi/2,'g',0,...
                      'D_A',0.,'rho_m',2,'nu',10,'Ca',0.02,...                       
                      'nParticles',0);
                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

