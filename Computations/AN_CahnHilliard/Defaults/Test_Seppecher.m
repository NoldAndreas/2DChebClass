                 
function [optsNum,optsPhys] = Test_Seppecher()

    %Numerical Parameters    
    Phys_Area = struct('y1Max',10,'N1',30,'N2',40); %y1Max = 20
    
    Plot_Area = struct('y1Min',0,'y1Max',Phys_Area.y1Max,'N1',100,...
                       'y2Min',0,'y2Max',pi,'N2',100);
                   
    %Sub_Area  = Phys_Area;
    Sub_Area = Phys_Area;
        
    optsNum = struct('PhysArea',Phys_Area,...
                     'PlotArea',Plot_Area,'SubArea',Sub_Area,...
                     'DDFTCode','Seppecher',...
                     'plotTimes',0:0.1:5);                            
                      
    optsPhys = struct('theta',pi/2,'g',0,...
                      'D_A',0,'D_B',0.,...
                      'rho_m',2,'nu',10,'Ca',0.02,... 
                     'kBT',0.7,...
                     'nParticles',0);
                 
    %Run File
    AddPaths();
    f = str2func(optsNum.DDFTCode);
    f(optsPhys,optsNum);                 

end                 

