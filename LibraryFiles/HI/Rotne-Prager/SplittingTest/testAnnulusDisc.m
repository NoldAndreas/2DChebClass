Phys_Area = struct('shape','InfSpace_FMT','y1Min',-inf,'y1Max',inf,'N',[20,20],'L1',10,...
                    'y2Min',-inf,'y2Max',inf,'L2',10);
                
IDC = InfSpace_FMT(Phys_Area);

IDC.ComputeConvolutionMatrix_Test_Short_3