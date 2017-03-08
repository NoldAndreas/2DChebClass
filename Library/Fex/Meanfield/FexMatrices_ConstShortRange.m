function convStruct = FexMatrices_ConstShortRange(optsPhys,IDC)

    lambda  = 1.5;
    LJsigma = 1;  
    params  = optsPhys.V2;
        
    shape.N      = optsPhys.FexNum.N;
    shape.volume = true;        
    shape.R      = LJsigma;

    conv1          = IDC.ComputeConvolutionFiniteSupport(Sphere(shape),'',IDC.Pts);
    
    shape.R = lambda;
    conv2 = IDC.ComputeConvolutionFiniteSupport(Sphere(shape),'',IDC.Pts);

    convStruct.Conv =  -(conv2(:,:,1) - conv1(:,:,1));
    
    %Rescale
    alpha           = -2/3*pi*(lambda^3 - LJsigma^3);        
    c               = (-16/9*pi)/alpha;    
    convStruct.Conv = c*convStruct.Conv;
    
    if(isfield(params,'epsilon'))
        convStruct.Conv = convStruct.Conv*params.epsilon;
    end

end