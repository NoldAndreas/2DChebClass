function FMTMatrices1DFull = Get1DMatrices(FMTMatrices3D,IDC,maskFull,maskComposed)
%function FMTMatrices = Polar_SpectralFourier_FMTMatrices(R,Pts,Maps,optsNum ,FFTMatrix)
    
    if(nargin == 2)
        maskFull     = (IDC.Pts.y1_kv    == inf);
        maskComposed = (IDC.AD.Pts.y1_kv == inf);            
    end
    
    for iSpecies = 1:length(IDC.R)
         
        
        FMTMatrices1D.AD.n3  = FMTMatrices3D.AD.n3(maskComposed,maskFull);
        FMTMatrices1D.AAD.n3 = FMTMatrices3D.AAD.n3(maskFull,maskComposed);

        FMTMatrices1D.AD.n2  = FMTMatrices3D.AD.n2(maskComposed,maskFull);
        FMTMatrices1D.AAD.n2 = FMTMatrices3D.AAD.n2(maskFull,maskComposed);
        
        %FMTMatrices1D.AD.n1  = FMTMatrices3D.AD.n1(maskComposed,maskFull);
        %FMTMatrices1D.AAD.n1 = FMTMatrices3D.AAD.n1(maskFull,maskComposed);

%        FMTMatrices1D.AD.n0  = FMTMatrices3D.AD.n0(maskComposed,maskFull);
%        FMTMatrices1D.AAD.n0 = FMTMatrices3D.AAD.n0(maskFull,maskComposed);

        %FMTMatrices1D.AD.n1_v_1  = FMTMatrices3D.AD.n1_v_1(maskComposed,maskFull);
        %FMTMatrices1D.AAD.n1_v_1 = FMTMatrices3D.AAD.n1_v_1(maskFull,maskComposed);
        
        %FMTMatrices1D.AD.n1_v_2  = FMTMatrices3D.AD.n1_v_2(maskComposed,maskFull);
        %FMTMatrices1D.AAD.n1_v_2 = FMTMatrices3D.AAD.n1_v_2(maskFull,maskComposed);

        FMTMatrices1D.AD.n2_v_1  = FMTMatrices3D.AD.n2_v_1(maskComposed,maskFull);
        FMTMatrices1D.AAD.n2_v_1 = FMTMatrices3D.AAD.n2_v_1(maskFull,maskComposed);

        FMTMatrices1D.AD.n2_v_2  = FMTMatrices3D.AD.n2_v_2(maskComposed,maskFull);
        FMTMatrices1D.AAD.n2_v_2 = FMTMatrices3D.AAD.n2_v_2(maskFull,maskComposed);

        
        FMTMatrices1DFull(iSpecies) = FMTMatrices1D;
    end
       
end

