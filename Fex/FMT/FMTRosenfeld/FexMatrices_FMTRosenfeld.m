function FMTMatricesFull = FexMatrices_FMTRosenfeld(Parameters,IDC)
%function FMTMatrices = Polar_SpectralFourier_FMTMatrices(R,Pts,Maps,optsNum ,FFTMatrix)
    
    v2struct(Parameters); 
    v2struct(FexNum);
    if(exist('sigmaS','var'))
        Rs  = diag(sigmaS)/2;    
    else
        Rs  = diag(V2.sigmaS)/2;    
    end
    
    %****************************
    %******* Test ********
    %****************************
    for iSpecies = 1:length(Rs)                           
        R        = Rs(iSpecies);
    
        shape.R  = R;
        shape.N  = [N1disc,N2disc];
                
        [AD,AAD] = IDC.GetAverageDensities(Disc(shape),{});

        FMTMatrices.AD.n2  = AD;
        FMTMatrices.AAD.n2 = AAD;    

        shape.N  = Ncircle;
        [AD,AAD] = IDC.GetAverageDensities(Circle(shape),{'cosPlc';'sinPlc'});

        FMTMatrices.AD.n1     = AD(:,:,1);
        FMTMatrices.AD.n1_v_1 = AD(:,:,2);
        FMTMatrices.AD.n1_v_2 = AD(:,:,3);

        FMTMatrices.AAD.n1     = AAD(:,:,1);
        FMTMatrices.AAD.n1_v_1 = AAD(:,:,2);
        FMTMatrices.AAD.n1_v_2 = AAD(:,:,3);

        FMTMatrices.AD.n0      = FMTMatrices.AD.n1/(2*pi*R);
        FMTMatrices.AAD.n0     = FMTMatrices.AAD.n1/(2*pi*R);                   
        
        FMTMatricesFull(iSpecies) = FMTMatrices;
    end
       
end

