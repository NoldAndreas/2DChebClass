function FMTMatricesFull = FexMatrices_FMTRosenfeld_3DFluid(Parameters,IDC)
%function FMTMatrices = Polar_SpectralFourier_FMTMatrices(R,Pts,Maps,optsNum ,FFTMatrix)
    global R 
    
    v2struct(Parameters);
    if(isfield(Parameters,'V2'))
        v2struct(V2);
    end
    v2struct(FexNum);
    
    Rs  = diag(sigmaS)/2;
    
    shape.N  = [N1disc,N2disc];	
    
    for iSpecies = 1:length(Rs)
         
        R        = Rs(iSpecies);
        shape.R  = R;
        
        shapeBall          = shape;
        shapeBall.theta1   = 0;
        shapeBall.theta2   = pi;                             

        %AAD                = IDC.AD.GetAverageDensities(Ball(shapeBall),{[-1/R 0],[0 -1/R]},IDC.Pts);
        [AD,AAD]           = IDC.GetAverageDensities(Ball(shapeBall),...
            {'weight_FMT_Roselfeld_x','weight_FMT_Roselfeld_y'});%{[-1/R 0],[0 -1/R]});        
        
        %shapeBall.N        = GetPointNumber(shapeBall,'Ball',requiredAbsoluteAccuracy);                
        %shapeBall.N        = [22,10];                
        
        FMTMatrices.AD.n2  = AD(:,:,1);
        FMTMatrices.AAD.n2 = AAD(:,:,1);    

        %FMTMatrices.AD.n0  = AD(:,:,1)/(4*pi*R^2);
        %FMTMatrices.AAD.n0 = AAD(:,:,1)/(4*pi*R^2);

        %FMTMatrices.AD.n1  = AD(:,:,1)/(4*pi*R);
        %FMTMatrices.AAD.n1 = AAD(:,:,1)/(4*pi*R);    

        FMTMatrices.AD.n2_v_1  = AD(:,:,2);
        FMTMatrices.AAD.n2_v_1 = AAD(:,:,2); 

        FMTMatrices.AD.n2_v_2 = AD(:,:,3);
        FMTMatrices.AAD.n2_v_2 = AAD(:,:,3); 

        %FMTMatrices.AD.n1_v_1  = AD(:,:,2)/(4*pi*R);
        %FMTMatrices.AAD.n1_v_1 = AAD(:,:,2)/(4*pi*R); 

        %FMTMatrices.AD.n1_v_2  = AD(:,:,3)/(4*pi*R);
        %FMTMatrices.AAD.n1_v_2 = AAD(:,:,3)/(4*pi*R);        
        
        
        shapeSphere         = shape; 
        shapeSphere.sphere  = true;
        
        %AAD                 = IDC.AD.GetAverageDensities(Disc(shapeSphere),'',IDC.Pts);
        %shapeSphere.N       = GetPointNumber(shapeSphere,'Disc',requiredAbsoluteAccuracy);        
        [AD,AAD]            = IDC.GetAverageDensities(Disc(shapeSphere),'');        

        FMTMatrices.AD.n3  = AD;
        FMTMatrices.AAD.n3 = AAD;
        
        FMTMatricesFull(iSpecies) = FMTMatrices;
    end
    

    

       
end

