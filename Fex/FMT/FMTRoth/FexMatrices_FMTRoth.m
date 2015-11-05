function FMTMatricesFull = FexMatrices_FMTRoth(Parameters,IDC)
    
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
        FMTMatrices.AAD.n2  = AAD;        

        shape.N  = Ncircle;
        
        [AD,AAD] = IDC.GetAverageDensities(Circle(shape), ...
                       {'cosPlc';'sinPlc';'cosSq';'sinSq';'cossin'});
                   
        FMTMatrices.AD.m0       = AD(:,:,1);
                   
        FMTMatrices.AD.m1_1     = AD(:,:,2);
        FMTMatrices.AD.m1_2     = AD(:,:,3);
        
        FMTMatrices.AD.m2_11    = AD(:,:,4);
        FMTMatrices.AD.m2_22    = AD(:,:,5);
        FMTMatrices.AD.m2_12    = AD(:,:,6);
        FMTMatrices.AD.m2_21    = FMTMatrices.AD.m2_12;
        
        FMTMatrices.AD.n0       = FMTMatrices.AD.m0/(2*pi*R);

                
        FMTMatrices.AAD.m0       = AAD(:,:,1);
                   
        FMTMatrices.AAD.m1_1     = AAD(:,:,2);
        FMTMatrices.AAD.m1_2     = AAD(:,:,3);
        
        FMTMatrices.AAD.m2_11    = AAD(:,:,4);
        FMTMatrices.AAD.m2_22    = AAD(:,:,5);
        FMTMatrices.AAD.m2_12    = AAD(:,:,6);
        FMTMatrices.AAD.m2_21    = FMTMatrices.AAD.m2_12;
        
        FMTMatrices.AAD.n0       = FMTMatrices.AAD.m0/(2*pi*R);
        
        FMTMatricesFull(iSpecies) = FMTMatrices;
        
    end    
             
end

